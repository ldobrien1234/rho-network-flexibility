# =====================================================================
#  cdc42_rac_rho.jl
#
#  A reaction-diffusion model of three Rho GTPases: Cdc42 (C), Rac (R)
#  and Rho (rho).
#
#  Kinetics.  Each GTPase is produced at a basal rate plus Hill-function
#  regulation by the others and by itself (autocatalysis), and decays
#  linearly:
#
#    dC/dt   = d_C   ΔC   + b_C   + a_C    H-(rho; K_C) + a_CC     Hs(C;   K_CC)     - g_C   C
#    dR/dt   = d_R   ΔR   + b_R   + a_R    H+(C;   K_R) + a_RR     Hs(R;   K_RR)     - g_R   R
#    drho/dt = d_rho Δrho + b_rho + a_rho1 H-(C;   K_rho1) + a_rho2 H+(R; K_rho2) + a_rhorho Hs(rho; K_rhorho) - g_rho rho
#
#  with Hill functions of order n (regulation) and n_s (self-activation)
#
#    H+(u; K) = u^n   / (K^n   + u^n)      (activation)
#    H-(u; K) = K^n   / (K^n   + u^n)      (inhibition)
#    Hs(u; K) = u^n_s / (K^n_s + u^n_s)    (self-activation)
#
#  So Cdc42 and Rho inhibit each other, Cdc42 activates Rac, and Rac
#  activates Rho.  The rates a, b and g are set per point by `lift` in
#  cycle_atlas.jl.
#
#  Bounded autocatalysis.  The self-activation terms a_ii Hs(u_i; K_ii)
#  let the diagonal of the Jacobian, f_ii = a_ii Hs'(u_i*) - g_i, be
#  positive, which the Hopf bands need.  Because they saturate and every
#  decay rate g_i is positive, all production is bounded while decay grows
#  linearly, so every solution stays bounded.  This enables oscillatory
#  solutions without blow-up.  (Getting f_ii > 0 from a negative decay rate
#  instead, as before, made the Hopf and mixed regions blow up.)
#
#  Positive concentrations.  No concentration can become negative if, at
#  u_i = 0, production is never negative, whatever the other species do:
#
#    b_i + sum over j != i of min(a_ij, 0) >= 0,
#
#  since each Hill function lies in [0, 1] and the self term is 0 at
#  u_i = 0.  (Diffusion with no-flux boundaries keeps u >= 0 by itself.)
#  `lift` chooses the decay rates g_i so that this holds.  It needs the
#  self-activation to be cooperative (n_s >= 2): with f_ii > 0 the self
#  term must be steeper than proportional at u*, Hs'(1) > Hs(1), and
#  u Hs'(u) / Hs(u) < n_s.
# =====================================================================

using NonlinearSolve
using DifferentialEquations

# --------------------------- model constants -------------------------
# The diagonal of the Jacobian, f_ii, is set in cycle_atlas.jl; the lift
# there chooses the self-activation amplitudes a_ii to produce it.
# With d = (1.8, 0.1, 0.15) and f = (-4.0, 1.4, 0.5) there are Hopf bands,
# but only with C_2 < 0.
# Baseline values with C_2 > 0 and steady-state bands only:
#   d = (0.05, 0.5, 0.05), f = (-5.0, -1.0, -5.0)

# Diffusion coefficients
const d_C   = 1.8
const d_R   = 0.1
const d_rho = 0.15

# The linear decay rates g_i are not constants: `lift` in cycle_atlas.jl
# chooses them per point, positive and large enough that concentrations
# stay positive (see the header).

# Hill order (n in the equations above).
# Chosen for stability of the PDE runs, not from the biology.  Tested
# 2026-11-25 with n = 1..4, when f_ii > 0 still came from negative decay
# rates: no n removed the blow-up, larger n blew up sooner, and n = 1 was
# the only value for which any Hopf-region run stayed bounded.
const n_hill = 1

# Hill thresholds
const K_C    = 0.5
const K_R    = 2.0
const K_rho1 = 0.5
const K_rho2 = 2.0

# Hill order of the self-activation (autocatalysis) terms (n_s in the
# equations above).  It must be at least 2 for concentrations to stay
# positive (see the header); 2 is the least cooperative choice.
const n_self = 2

# Hill thresholds of the self-activation terms
const K_CC     = 2.0
const K_RR     = 2.0
const K_rhorho = 2.0

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

"Self-activation Hill function Hs(u; K) = u^n_s / (K^n_s + u^n_s)."
function hill_self(u, K)
    u = max(u, 0.0)
    return u^n_self / (K^n_self + u^n_self)
end

"Derivative of the self-activation Hill function, dHs/du = n_s K^n_s u^(n_s-1) / (K^n_s + u^n_s)^2."
function dhill_self(u, K)
    u = max(u, 0.0)
    return n_self * K^n_self * u^(n_self - 1) / (K^n_self + u^n_self)^2
end

# --------------------------- reaction term ---------------------------
"""
    reaction(u, p)

Reaction term F(u) for u = (C, R, rho), the right-hand sides of the
equations in the header.

`p` holds the rates that are not fixed constants:
    p.a_C, p.a_R, p.a_rho1, p.a_rho2    Hill amplitudes of the regulation
    p.a_CC, p.a_RR, p.a_rhorho          Hill amplitudes of the self-activation
    p.b_C, p.b_R, p.b_rho               basal production rates
    p.g_C, p.g_R, p.g_rho               linear decay rates
"""
function reaction(u, p)
    C, R, rho = u

    # regulation by the other GTPases
    regulation_C   = p.a_C    * hill_inhibit(rho, K_C)
    regulation_R   = p.a_R    * hill_activate(C, K_R)
    regulation_rho = p.a_rho1 * hill_inhibit(C, K_rho1) + p.a_rho2 * hill_activate(R, K_rho2)

    # bounded autocatalysis: saturating, so it cannot outgrow the linear decay
    self_C   = p.a_CC     * hill_self(C,   K_CC)
    self_R   = p.a_RR     * hill_self(R,   K_RR)
    self_rho = p.a_rhorho * hill_self(rho, K_rhorho)

    F_C   = p.b_C   + regulation_C   + self_C   - p.g_C   * C
    F_R   = p.b_R   + regulation_R   + self_R   - p.g_R   * R
    F_rho = p.b_rho + regulation_rho + self_rho - p.g_rho * rho

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
# the resulting ODEs are integrated in time (`simulate`).

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
    simulate(u_init, p; t_max = 2_000.0, tol = 1e-8, blow_up = 1e3, save_every = 5.0)

Integrate the PDE in time from the initial value `u_init` (an N x 3 matrix
with columns C, R, rho on N equally spaced points of [0, L]).  The same
integration is used for every point, whatever the linear theory predicts,
and the outcome is classified afterwards.  It stops at the first of:

    :steady_state   the solution stops changing (|du/dt| below `tol`)
    :blow_up        some concentration exceeds `blow_up` in size
    :not_settled    time `t_max` is reached: the solution may oscillate,
                    or still be on its way to a steady state
    :failed         the solver stops for any other reason

With every decay rate g_i > 0 (bounded autocatalysis, see the header)
solutions stay bounded, so :blow_up should not occur; it is kept as a
safeguard for decay rates g_i < 0, where nothing balances the growth
-g_i u_i once the Hill terms saturate.

`tol` is also the solver's absolute and relative tolerance.  It must be
small: patterns can sit near an unstable state for a long time before
moving on, and a loose tolerance can stop there and report it as the
steady state.

Returns a named tuple
    status      one of the four outcomes above
    times       the times of the saved snapshots: every `save_every`,
                and the time the integration stopped
    snapshots   the N x 3 solution at each of `times`
"""
function simulate(u_init, p; t_max = 2_000.0, tol = 1e-8, blow_up = 1e3, save_every = 5.0)
    # stop when the solution stops changing
    steady_state_callback = TerminateSteadyState(tol, tol)

    # stop when some concentration exceeds blow_up
    is_blowing_up(u, t, integrator) = maximum(abs, u) > blow_up
    blow_up_callback = DiscreteCallback(is_blowing_up, terminate!)

    problem  = ODEProblem(pde_rhs, u_init, (0.0, t_max), p)
    solution = solve(problem, FBDF(); abstol = tol, reltol = tol, saveat = save_every,
                     callback = CallbackSet(steady_state_callback, blow_up_callback))

    # both callbacks end with the return code Terminated; the final state tells them apart
    u_end = solution.u[end]
    if solution.retcode == ReturnCode.Terminated
        status = maximum(abs, u_end) > blow_up ? :blow_up : :steady_state
    elseif solution.retcode == ReturnCode.Success
        status = :not_settled
    else
        status = :failed
    end

    return (status = status, times = solution.t, snapshots = solution.u)
end




