# =====================================================================
#  cdc42_rac_rho.jl
#
#  A reaction-diffusion model of three Rho GTPases: Cdc42 (C), Rac (R)
#  and Rho (rho).
#
#  Kinetics.  Each GTPase is produced at a basal rate plus Hill-function
#  regulation by the others, and decays linearly:
#
#    dC/dt   = d_C * ΔC + b_C   + a_C   H-(rho; K_C) - g_C   C
#    dR/dt   = d_R * ΔR + b_R   + a_R   H+(C;   K_R) - g_R   R
#    drho/dt =d_rho * Δrho +  b_rho + a_rho1 H-(C;  K_rho1) + a_rho2 H+(R; K_rho2) - g_rho rho
#
#  with Hill functions of order n
#
#    H+(u; K) = u^n / (K^n + u^n)     (activation)
#    H-(u; K) = K^n / (K^n + u^n)     (inhibition)
#
#  So Cdc42 and Rho inhibit each other, Cdc42 activates Rac, and Rac
#  activates Rho.
# =====================================================================

using NonlinearSolve
using DifferentialEquations

# --------------------------- model constants -------------------------
# Values from the original cycle_atlas.jl.  They give steady-state bands
# with C_2 > 0 (mutual inhibition of Cdc42 and Rho) and no Hopf bands.
# Values that do give Hopf bands, but only with C_2 < 0, are
#   d = (1.8, 0.1, 0.15),  g = (4.0, -1.4, -0.5)
#Baseline values with C_2 > 0
#   d = (0.05, 0.5, 0.05), g = (5.0, 1.0, 5.0)
# (negative decay rates for R and rho: net self-activation).

# Diffusion coefficients
const d_C   = 1.8
const d_R   = 0.1
const d_rho = 0.15

# Linear decay rates
const g_C   = 4.0
const g_R   = -1.4
const g_rho = -0.5

# Hill order (n in the equations above)
const n_hill = 1

# Hill thresholds
const K_C    = 0.5
const K_R    = 2.0
const K_rho1 = 0.5
const K_rho2 = 2.0

# Domain [0, L] with no-flux boundaries
const L = π

# --------------------------- Hill functions --------------------------
# Negative concentrations are treated as zero, so that u^n stays >= 0
# and K^n + u^n cannot vanish.

"Activating Hill function H+(u; K) = u^n / (K^n + u^n)."
function hill_activate(u, K)
    u = max(u, 0.0) #must take positive u values
    return u^n_hill / (K^n_hill + u^n_hill)
end

"Inhibiting Hill function H-(u; K) = K^n / (K^n + u^n)."
function hill_inhibit(u, K)
    u = max(u, 0.0) #must take positive u values
    return K^n_hill / (K^n_hill + u^n_hill)
end

"Derivative of the activating Hill function, dH+/du = n K^n u^(n-1) / (K^n + u^n)^2."
function dhill_activate(u, K)
    u = max(u, 0.0)
    return n_hill * K^n_hill * u^(n_hill - 1) / (K^n_hill + u^n_hill)^2
end

"Derivative of the inhibiting Hill function, dH-/du = -n K^n u^(n-1) / (K^n + u^n)^2."
function dhill_inhibit(u, K)
    u = max(u, 0.0)
    return -n_hill * K^n_hill * u^(n_hill - 1) / (K^n_hill + u^n_hill)^2
end

# --------------------------- reaction term ---------------------------
"""
    reaction(u, p)

Reaction term F(u) for u = (C, R, rho), the right-hand sides of the
equations in the header.

`p` holds the rates that are not fixed constants:
    p.a_C, p.a_R, p.a_rho1, p.a_rho2    Hill amplitudes
    p.b_C, p.b_R, p.b_rho               basal production rates
"""
function reaction(u, p)
    C, R, rho = u

    F_C   = p.b_C   + p.a_C    * hill_inhibit(rho, K_C)                                        - g_C   * C
    F_R   = p.b_R   + p.a_R    * hill_activate(C, K_R)                                         - g_R   * R
    F_rho = p.b_rho + p.a_rho1 * hill_inhibit(C, K_rho1) + p.a_rho2 * hill_activate(R, K_rho2) - g_rho * rho

    return [F_C, F_R, F_rho]
end

# --------------------------- homogeneous steady state ----------------
"""
    homogeneous_state(p; u0 = [1.0, 1.0, 1.0])

Solve F(u*) = 0 for a spatially uniform steady state u* = (C*, R*, rho*),
starting from `u0` using the root finder NonlinearProblem
"""
function homogeneous_state(p; u0 = [1.0, 1.0, 1.0])
    problem  = NonlinearProblem((u, p) -> reaction(u, p), u0, p)
    solution = solve(problem)
    u_star   = solution.u

    #Error messages if solver fails or finds a negative solution
    SciMLBase.successful_retcode(solution) ||
        error("no steady state found from u0 = $u0 (retcode $(solution.retcode))")
    all(u_star .> 0) ||
        error("steady state $u_star found from u0 = $u0 is not positive")

    return u_star
end

# --------------------------- PDE solver ------------------------------
# The PDE  du/dt = D u_xx + F(u)  on [0, L] with no-flux boundaries is
# discretised in space on N equally spaced points (method of lines), and
# the resulting ODEs are integrated in time until they reach steady state.

"""
    laplacian_neumann(v, h)

Second-order finite-difference approximation of v_xx on an equally spaced
grid with spacing h and no-flux boundaries.  At each end the ghost point
outside the domain is set equal to its mirror image inside (v_0 = v_2 and
v_{N+1} = v_{N-1}), which makes v_x = 0 there.
"""
function laplacian_neumann(v, h)
    #calculates second spatial derivative with finite difference approximation
    N = length(v)
    v_xx = similar(v)
    v_xx[1] = 2 * (v[2] - v[1]) / h^2
    for i in 2:N-1
        v_xx[i] = (v[i-1] - 2 * v[i] + v[i+1]) / h^2
    end
    v_xx[N] = 2 * (v[N-1] - v[N]) / h^2
    return v_xx
end

"""
    pde_rhs(u, p, t)

Right-hand side D u_xx + F(u) of the discretised PDE.  `u` is an N x 3
matrix whose columns are C, R and rho on the grid.
"""
function pde_rhs(u, p, t)
    N = size(u, 1)
    h = L / (N - 1)

    #diffusion term
    du = zeros(eltype(u), N, 3)
    du[:, 1] = d_C   * laplacian_neumann(u[:, 1], h)
    du[:, 2] = d_R   * laplacian_neumann(u[:, 2], h)
    du[:, 3] = d_rho * laplacian_neumann(u[:, 3], h)

    #reaction term
    for i in 1:N
        du[i, :] += reaction(u[i, :], p)
    end
    return du
end

"""
    solve_to_steady_state(u_init, p; t_max = 10_000.0, tol = 1e-8)

Integrate the PDE from the initial value `u_init` (an N x 3 matrix with
columns C, R, rho on N equally spaced points of [0, L]) until it stops
changing, using DifferentialEquations' SteadyStateProblem.  Gives up at
time `t_max`.

`tol` is used as both the absolute and relative tolerance.  It must be
small: patterns can sit near an unstable state for a long time before
moving on, and a loose tolerance (such as the solver's default) can stop
there and report it as the steady state.

Returns the final N x 3 matrix and whether a steady state was reached.
"""
function solve_to_steady_state(u_init, p; t_max = 10_000.0, tol = 1e-8)
    problem  = SteadyStateProblem(pde_rhs, u_init, p)
    solution = solve(problem, DynamicSS(FBDF(); tspan = t_max); abstol = tol, reltol = tol)

    reached = SciMLBase.successful_retcode(solution)
    return solution.u, reached
end




