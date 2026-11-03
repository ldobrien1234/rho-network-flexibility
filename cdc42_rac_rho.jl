# =====================================================================
#  cdc42_rac_rho.jl
#
#  A standard three-node reaction-diffusion model of the reduced
#  Cdc42 (C) -- Rac (R) -- Rho (rho) network analysed in the companion
#  LaTeX section "Flexibility of Spatial Patterning in the
#  Cdc42--Rac--Rho Network".
#
#  Interaction graph (reduced network, direct rho->R edge omitted):
#
#        C  --(+)-->  R  --(s)-->  rho
#        ^                          |
#        |                          |
#        +--------(-) / (-)---------+     (C <-> rho mutual inhibition)
#
#  with s = + in the PATTERNING regime (P): Rac activates Rho   (f_rhoR > 0)
#       s = - in the SWITCH    regime (S): Rac inhibits Rho     (f_rhoR < 0)
#
#  Literature basis for the two signs of the Rac -> Rho edge:
#    * Rac INHIBITS Rho: active Rac -> Pak1/Pak4 -> phosphorylation and
#      inactivation of RhoA-directed GEFs (p115-RhoGEF, GEF-H1, PDZ-RhoGEF).
#      Sanz-Moreno et al.; Rosenfeldt et al. 2006 (J Mol Signal 1:8);
#      Guilluy et al. 2011 (Nat Cell Biol review "Rho protein crosstalk").
#    * Rac ACTIVATES Rho: active Rac recruits the Lbc-type GEFs
#      Arhgef11/Arhgef12 as effectors, driving Rho during retraction.
#      Marston et al. 2023, Nat Commun 14:8151
#      ("Rho GTPase activity crosstalk mediated by Arhgef11 and Arhgef12").
#    * Switch / winner-take-all polarity (regime S phenomenology):
#      Mori, Jilkine & Edelstein-Keshet 2008, Biophys J 94:3684 (wave-pinning);
#      Jilkine & Edelstein-Keshet 2011, PLoS Comput Biol (GTPase polarity review).
#
#  Kinetics: each node is a single activity variable with linear turnover
#  and additive Hill regulation, the minimal smooth realisation consistent
#  with the interaction graph and with f_ii < 0.  Written so that the
#  reaction Jacobian at a prescribed homogeneous state can be set exactly,
#  which is what makes the constructive Proposition of the text usable.
# =====================================================================

using LinearAlgebra
using ForwardDiff
using Random
using Printf

# --------------------------- Hill functions ---------------------------
# activating (increasing 0..1) and inhibiting (decreasing 1..0) Hill terms
@inline hillA(u, K, m) = (uc = max(u, 0.0); um = uc^m; km = K^m; um / (km + um))
@inline hillI(u, K, m) = (uc = max(u, 0.0); um = uc^m; km = K^m; km / (km + um))
@inline dhillA(u, K, m) = (uc = max(u, 1e-12); m * K^m * uc^(m - 1) / (K^m + uc^m)^2)
@inline dhillI(u, K, m) = -dhillA(u, K, m)

# --------------------------- parameter container ---------------------
"""
Parameters of the three-node RD model.

Reaction terms (u = (C,R,rho)):
  f_C   = bC  + aC  * hillI(rho, KC,   m)                          - gC  * C
  f_R   = bR  + aR  * hillA(C,   KR,   m)                          - gR  * R
  f_rho = brho+ arho1*hillI(C,   Krho1,m) + srho*arho2*hillX(R,Krho2,m) - grho* rho

  regime == :P  ->  hillX = hillA, srho = +1   (Rac activates Rho)
  regime == :S  ->  hillX = hillI, srho = +1   (Rac inhibits Rho; term is decreasing)
"""
Base.@kwdef struct RDParams
    dC::Float64;  dR::Float64;  drho::Float64
    gC::Float64;  gR::Float64;  grho::Float64
    m::Int
    bC::Float64;  aC::Float64;  KC::Float64
    bR::Float64;  aR::Float64;  KR::Float64
    brho::Float64; arho1::Float64; Krho1::Float64
    arho2::Float64; Krho2::Float64
    regime::Symbol            # :P or :S
    label::String = ""
    target_mode::Int = -1     # design target (for bookkeeping)
end

diffusion(p::RDParams) = (p.dC, p.dR, p.drho)

# --------------------------- reaction kinetics -----------------------
"local reaction term F(u); u is a 3-vector (C,R,rho)."
function reaction(u, p::RDParams)
    C, R, rho = u
    hillX = p.regime === :P ? hillA : hillI
    fC   = p.bC   + p.aC    * hillI(rho, p.KC,   p.m)                              - p.gC   * C
    fR   = p.bR   + p.aR    * hillA(C,   p.KR,   p.m)                              - p.gR   * R
    frho = p.brho + p.arho1 * hillI(C,   p.Krho1, p.m) + p.arho2 * hillX(R, p.Krho2, p.m) - p.grho * rho
    return SVectorlike(fC, fR, frho)
end

# lightweight fixed 3-vector to avoid a StaticArrays dependency
SVectorlike(a, b, c) = [a, b, c]

"reaction Jacobian at u (3x3), via automatic differentiation."
reaction_jac(u, p::RDParams) = ForwardDiff.jacobian(v -> reaction(v, p), collect(float.(u)))

# --------------------------- homogeneous steady state ---------------
"""
Newton solve for a homogeneous steady state near `u0`.
Returns (u_star, converged::Bool).
"""
function homogeneous_state(p::RDParams; u0 = [1.0, 1.0, 1.0], tol = 1e-11, maxit = 100)
    u = collect(float.(u0))
    for _ in 1:maxit
        F = reaction(u, p)
        nrm = norm(F)
        nrm < tol && return (u, true)
        J = reaction_jac(u, p)
        du = J \ F
        # damped Newton, keep positivity
        α = 1.0
        while α > 1e-6
            un = u .- α .* du
            if all(un .> 0) && norm(reaction(un, p)) < nrm
                u = un; break
            end
            α *= 0.5
        end
        α <= 1e-6 && (u = max.(u .- du, 1e-9))
    end
    return (u, norm(reaction(u, p)) < 1e-7)
end

# --------------------------- linear dispersion ----------------------
"""
Growth rate of admissible mode n (wavenumber kappa = n on [0,pi], B = -n^2):
returns (max_real_part, is_complex_pair::Bool, detL::Float64) for L(n) = J - n^2 D.
"""
function mode_spectrum(n::Integer, J::AbstractMatrix, D::NTuple{3,Float64})
    L = J - (n^2) .* Diagonal(collect(D))
    ev = eigvals(L)
    i = argmax(real.(ev))
    λ = ev[i]
    iscplx = abs(imag(λ)) > 1e-9
    return (real(λ), iscplx, real(det(L)))
end

"""
Dispersion data over modes 0:nmax at a given homogeneous state.
Returns a NamedTuple with vectors `n`, `growth`, `complex`, `detL`,
plus `homog_stable`, the fastest-growing mode `nstar` with its `growth_star`,
and `turing` (true iff some n>=1 has a positive real eigenvalue).
"""
function dispersion(p::RDParams; nmax = 14, ustar = nothing)
    us = ustar === nothing ? homogeneous_state(p)[1] : ustar
    J  = reaction_jac(us, p)
    D  = diffusion(p)
    ns = collect(0:nmax)
    growth  = similar(ns, Float64)
    cplx    = falses(length(ns))
    detLv   = similar(ns, Float64)
    for (k, n) in enumerate(ns)
        g, c, d = mode_spectrum(n, J, D)
        growth[k] = g; cplx[k] = c; detLv[k] = d
    end
    homog_stable = growth[1] < 1e-9          # n = 0
    # fastest growing among n >= 1
    idx = 2:length(ns)
    j   = idx[argmax(growth[idx])]
    any_pos_mode = any((growth[k] > 1e-7) && !cplx[k] for k in idx)
    # a *finite-wavelength* (Turing) instability: homogeneous stable, some n>=1 unstable
    turing = homog_stable && any_pos_mode
    return (; n = ns, growth, complex = cplx, detL = detLv,
              growth0 = growth[1], homog_stable, nstar = ns[j], growth_star = growth[j],
              any_pos_mode, turing, ustar = us, J)
end

# =====================================================================
#  1-D PDE solver:  u_t = D u_xx + F(u)  on [0,pi], no-flux (modes cos(n x), B = -n^2)
#  Method of lines + IMEX Euler (implicit diffusion, explicit reaction).
# =====================================================================
"""
    simulate(p; N, L, tmax, dt, seed, amp, tol)

Integrate to (numerical) steady state from a small random perturbation of
the homogeneous state.  Returns a NamedTuple with the grid `x`, the final
profiles `C,R,rho`, the homogeneous state `ustar`, the stopping time `t`,
and the final max |du/dt| `resid`.
"""
function simulate(p::RDParams; N = 257, L = π, tmax = 4000.0, dt = 4e-3,
                  seed = 1, amp = 1e-2, tol = 1e-8, ic = :perturb, seed_mode = 0)
    x  = range(0, L; length = N)
    h  = step(x)
    us, ok = homogeneous_state(p)
    ok || @warn "homogeneous_state did not converge for $(p.label)"

    # Neumann Laplacian (second order), as a Tridiagonal
    dl = fill(1.0 / h^2, N - 1)
    d0 = fill(-2.0 / h^2, N)
    du = fill(1.0 / h^2, N - 1)
    du[1]      = 2.0 / h^2          # ghost node u_{-1} = u_{1}
    dl[end]    = 2.0 / h^2          # ghost node u_{N}  = u_{N-2}
    Lap = Tridiagonal(dl, d0, du)
    I_N = Diagonal(ones(N))
    facs = (lu(I_N - dt * p.dC   * Lap),
            lu(I_N - dt * p.dR   * Lap),
            lu(I_N - dt * p.drho * Lap))

    rng = _rng(seed)
    C   = fill(us[1], N); R = fill(us[2], N); rho = fill(us[3], N)
    if ic === :perturb
        C   .+= amp .* us[1] .* randn(rng, N)
        R   .+= amp .* us[2] .* randn(rng, N)
        rho .+= amp .* us[3] .* randn(rng, N)
    elseif ic === :mode            # seed the prescribed admissible mode + tiny noise
        cx = cos.(seed_mode .* x)
        C   .+= amp .* us[1] .* cx .+ 1e-3 .* us[1] .* randn(rng, N)
        R   .+= amp .* us[2] .* cx .+ 1e-3 .* us[2] .* randn(rng, N)
        rho .-= amp .* us[3] .* cx .+ 1e-3 .* us[3] .* randn(rng, N)
    elseif ic === :step            # winner-take-all probe: half domain high
        half = N ÷ 2
        C[1:half]   .*= 1.6;  C[half+1:end]   .*= 0.5
        rho[1:half] .*= 0.5;  rho[half+1:end] .*= 1.6
    end
    C .= max.(C, 1e-9); R .= max.(R, 1e-9); rho .= max.(rho, 1e-9)

    nsteps = ceil(Int, tmax / dt)
    t, resid = _evolve!(C, R, rho, facs, p.regime === :P,
                        p.bC, p.aC, p.KC, p.gC,
                        p.bR, p.aR, p.KR, p.gR,
                        p.brho, p.arho1, p.Krho1, p.arho2, p.Krho2, p.grho,
                        p.m, dt, nsteps, tol)
    return (; x = collect(x), C, R, rho, ustar = us, t, resid,
              amp_final = maximum(rho) - minimum(rho))
end

# typed inner loop (function barrier): IMEX Euler steps to steady state
function _evolve!(C::Vector{Float64}, R::Vector{Float64}, rho::Vector{Float64},
                  facs, isP::Bool,
                  bC, aC, KC, gC, bR, aR, KR, gR,
                  brho, arho1, Krho1, arho2, Krho2, grho,
                  m::Int, dt::Float64, nsteps::Int, tol::Float64)
    N = length(C)
    rC = similar(C); rR = similar(R); rrho = similar(rho); tmp = similar(C)
    Cp = copy(C); Rp = copy(R); rhop = copy(rho)
    fC, fR, frho = facs
    t = 0.0; resid = Inf
    for k in 1:nsteps
        @inbounds @simd for i in 1:N
            gcr = arho1 * hillI(C[i], Krho1, m) +
                  arho2 * (isP ? hillA(R[i], Krho2, m) : hillI(R[i], Krho2, m))
            rC[i]   = bC   + aC * hillI(rho[i], KC, m) - gC   * C[i]
            rR[i]   = bR   + aR * hillA(C[i],   KR, m) - gR   * R[i]
            rrho[i] = brho + gcr                       - grho * rho[i]
        end
        @inbounds @simd for i in 1:N
            tmp[i] = C[i] + dt * rC[i]
        end
        ldiv!(fC, tmp); copyto!(C, tmp)
        @inbounds @simd for i in 1:N
            tmp[i] = R[i] + dt * rR[i]
        end
        ldiv!(fR, tmp); copyto!(R, tmp)
        @inbounds @simd for i in 1:N
            tmp[i] = rho[i] + dt * rrho[i]
        end
        ldiv!(frho, tmp); copyto!(rho, tmp)
        @inbounds @simd for i in 1:N
            C[i]   = ifelse(C[i]   < 1e-12, 1e-12, C[i])
            R[i]   = ifelse(R[i]   < 1e-12, 1e-12, R[i])
            rho[i] = ifelse(rho[i] < 1e-12, 1e-12, rho[i])
        end
        t += dt
        if k % 100 == 0
            m1 = 0.0
            @inbounds for i in 1:N
                m1 = max(m1, abs(C[i]-Cp[i]), abs(R[i]-Rp[i]), abs(rho[i]-rhop[i]))
            end
            resid = m1 / (100dt)          # max |du/dt| over last 100 steps
            copyto!(Cp, C); copyto!(Rp, R); copyto!(rhop, rho)
            resid < tol && break
        end
    end
    return t, resid
end

# reproducible RNG
_rng(seed) = Random.Xoshiro(seed)
