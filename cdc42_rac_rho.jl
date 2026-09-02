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
Growth rate of admissible mode n (wavenumber kappa = n on [0,2pi], B = -n^2):
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
#  Constructive parameter synthesis  (Proposition "Existence of a
#  primary instability at a prescribed mode" in the text)
# =====================================================================
"""
    construct_regime_P(ntarget; ...)

Build a regime-(P) parameter set whose dispersion cubic g(B) = det L(B)
has a DOUBLE root at B = A = -ntarget^2 (marginal / codimension-one), with
the third root B3 > 0 and every other admissible mode strictly stable, and
with a linearly stable homogeneous state.  Implements the backward
construction of the text: fix the desired root geometry of g, solve for the
cycle quantities C_{Crho} (2-cycle) and C_3 (3-cycle), then realise them
with sign-legal Hill kinetics.

A subsequent call to `push_to_onset` moves a single 2-cycle gain slightly
past the marginal point to obtain a genuine one-mode instability.
"""
function construct_regime_P(ntarget::Integer;
        dR = 10.0, g = 1.0, m = 3,
        ustar = (1.0, 1.0, 1.0),
        KC = 0.2, Krho1 = 0.2, KR = 5.0, Krho2 = 5.0,
        dc_ladder = (0.2, 0.1, 0.06, 0.04, 0.025, 0.015, 0.009,
                     0.005, 0.003, 0.0018, 0.001, 0.0006, 0.0004))

    A   = -float(ntarget)^2
    fCC = -g; fRR = -g; frr = -g          # prescribed self-decays (f_ii = -gamma)
    uC, uR, urho = ustar

    for dc in dc_ladder
        dC = dc; drho = dc
        a3 = dC * dR * drho
        c2 = dC*dR*frr + dC*drho*fRR + dR*drho*fCC
        B3 = -c2 / a3 - 2A
        B3 > 0 || continue
        c1s = a3 * A * (A + 2B3)                    # target c1  (< 0)
        c0s = -a3 * A^2 * B3                        # target c0  (< 0)
        P1  = dC*fRR*frr + dR*fCC*frr + drho*fCC*fRR
        P0  = fCC * fRR * frr
        CCrho = (P1 - c1s) / dR                     # 2-cycle  f_Crho f_rhoC  (> 0)
        C3    = c0s - P0 + fRR * CCrho              # 3-cycle  f_RC f_rhoR f_Crho  (< 0)
        (CCrho > 0 && C3 < 0) || continue

        # homogeneous (n = 0) Routh-Hurwitz check with these Jacobian entries
        trJ = fCC + fRR + frr
        pJ  = fCC*fRR + fCC*frr + fRR*frr - CCrho   # sum of 2x2 principal minors
        detJ = c0s
        (trJ < 0 && detJ < 0 && trJ*pJ < detJ && pJ > 0) || continue

        # ---- realise the two cycle products with sign-legal edges ----
        # free splitting parameters theta (|f_Crho|) and phi (f_RC)
        best = nothing
        for θ in 0.4:0.1:2.9, φ in 0.4:0.1:2.9
            fCrho = -θ
            frhoC = -CCrho / θ
            fRC   =  φ
            frhoR =  C3 / (fRC * fCrho)             # = -C3/(phi*theta) > 0
            aC   = fCrho / dhillI(urho, KC,    m)   # (<0)/(<0) > 0
            arho1 = frhoC / dhillI(uC,   Krho1, m)
            aR   = fRC   / dhillA(uC,   KR,    m)
            arho2 = frhoR / dhillA(uR,   Krho2, m)
            (aC > 0 && arho1 > 0 && aR > 0 && arho2 > 0) || continue
            bC   = g*uC   - aC   * hillI(urho, KC,    m)
            bR   = g*uR   - aR   * hillA(uC,   KR,    m)
            brho = g*urho - arho1* hillI(uC,   Krho1, m) - arho2 * hillA(uR, Krho2, m)
            mb = min(bC, bR, brho)
            if best === nothing || mb > best.mb
                best = (; θ, φ, aC, arho1, aR, arho2, bC, bR, brho, mb,
                          dC, drho, CCrho, C3, B3, c1s, c0s)
            end
        end
        best === nothing && continue
        best.mb > 0.02 || continue

        p = RDParams(; dC = best.dC, dR = dR, drho = best.drho,
              gC = g, gR = g, grho = g, m = m,
              bC = best.bC, aC = best.aC, KC = KC,
              bR = best.bR, aR = best.aR, KR = KR,
              brho = best.brho, arho1 = best.arho1, Krho1 = Krho1,
              arho2 = best.arho2, Krho2 = Krho2,
              regime = :P, target_mode = ntarget,
              label = "P: primary Turing mode n=$(ntarget)")

        info = (; CCrho = best.CCrho, C3 = best.C3, B3 = best.B3,
                  c1 = best.c1s, c0 = best.c0s, dC = best.dC)
        return (p, info)
    end
    error("construct_regime_P: no feasible synthesis for n = $ntarget " *
          "(try widening dc_ladder or relaxing K/ustar).")
end

"""
    apply_twocycle_gain(p, gain; ustar)

Scale the C<->rho 2-cycle strength by `gain` (multiply aC and arho1 by sqrt),
recomputing bC, brho so the homogeneous state stays at `ustar`.  `gain = 1`
returns the marginal set; `gain` slightly above 1 is a small supercriticality.
"""
function apply_twocycle_gain(p::RDParams, gain; ustar = (1.0, 1.0, 1.0))
    uC, uR, urho = ustar
    s = sqrt(gain)
    aC    = p.aC    * s
    arho1 = p.arho1 * s
    bC   = p.gC  *uC   - aC   * hillI(urho, p.KC,    p.m)
    brho = p.grho*urho - arho1* hillI(uC,   p.Krho1, p.m) - p.arho2*hillA(uR, p.Krho2, p.m)
    RDParams(; p.dC, p.dR, p.drho, p.gC, p.gR, p.grho, p.m,
               bC, aC, p.KC, p.bR, p.aR, p.KR,
               brho, arho1, p.Krho1, p.arho2, p.Krho2,
               regime = :P, target_mode = p.target_mode, label = p.label)
end

"""
    push_to_onset(p; gains, gmin)

Search a rising 2-cycle gain for the smallest supercriticality at which the
target mode is the unique unstable admissible mode with growth > `gmin`.
Returns (p, gain, dispersion).  (Kept for exploration; the run script uses a
fixed small gain via `apply_twocycle_gain`.)
"""
function push_to_onset(p::RDParams; ustar = (1.0, 1.0, 1.0),
                       gains = 1.0:0.0004:1.06, gmin = 2.5e-3)
    mk = gain -> apply_twocycle_gain(p, gain; ustar)
    fallback = nothing
    for gain in gains
        q = mk(gain)
        d = dispersion(q; nmax = 16)
        nun = count(>(1e-7), d.growth[2:end])
        if d.homog_stable && d.nstar == p.target_mode && d.growth_star > gmin
            fallback === nothing && (fallback = (q, gain, d))
            nun == 1 && return (q, gain, d)          # unique unstable admissible mode
        end
    end
    return fallback === nothing ? (mk(gains[end]), gains[end], dispersion(mk(gains[end]))) : fallback
end

"""
    flip_to_regime_S(p)

Return the SAME parameter set with the single Rac->Rho edge flipped from
activating to inhibiting (regime S), recomputing brho so the homogeneous
state is unchanged.  This is the "delete/negate one regulatory edge"
operation of Corollary (Regime dichotomy).
"""
function flip_to_regime_S(p::RDParams; ustar = (1.0, 1.0, 1.0), edge_scale = 1.0,
                          note = "Rac->Rho edge negated")
    uC, uR, urho = ustar
    arho2 = p.arho2 * edge_scale
    brho  = p.grho*urho - p.arho1*hillI(uC, p.Krho1, p.m) - arho2*hillI(uR, p.Krho2, p.m)
    RDParams(; p.dC, p.dR, p.drho, p.gC, p.gR, p.grho, p.m,
               p.bC, p.aC, p.KC, p.bR, p.aR, p.KR,
               brho, p.arho1, p.Krho1, arho2, p.Krho2,
               regime = :S, target_mode = p.target_mode,
               label = replace(p.label, "P:" => "S:") * "  ($note)")
end

"""
    build_regime_S(; strength, ...)

Directly build a regime-(S) parameter set (Rac inhibits Rho) with prescribed
homogeneous state u*=(1,1,1) and sign-legal Hill kinetics, with all four
regulatory edges scaled by `strength`.  Used to show that NO choice of
`strength` produces a finite-wavelength instability: small `strength` keeps
a stable homogeneous state (perturbations decay); large `strength` makes the
network bistable (a winner-take-all switch) -- never a Turing pattern.
"""
function build_regime_S(; strength = 1.0, dC = 0.1, dR = 10.0, drho = 0.1,
        g = 1.0, m = 3, KC = 1.0, Krho1 = 1.0, KR = 5.0, Krho2 = 1.0,
        q = 1.0, r = 1.0, w = 1.0, label = "S: Rac inhibits Rho")
    uC = uR = urho = 1.0
    fCrho = -q * strength                  # Rho -| Cdc42
    frhoC = -q * strength                  # Cdc42 -| Rho
    fRC   =  r * strength                  # Cdc42 -> Rac
    frhoR = -w * strength                  # Rac -| Rho   (regime S: negative)
    aC    = fCrho / dhillI(urho, KC,    m)
    arho1 = frhoC / dhillI(uC,   Krho1, m)
    aR    = fRC   / dhillA(uC,   KR,    m)
    arho2 = frhoR / dhillI(uR,   Krho2, m)      # dhillI < 0, frhoR < 0  -> arho2 > 0
    bC   = g*uC   - aC   * hillI(urho, KC,    m)
    bR   = g*uR   - aR   * hillA(uC,   KR,    m)
    brho = g*urho - arho1* hillI(uC,   Krho1, m) - arho2 * hillI(uR, Krho2, m)
    RDParams(; dC, dR, drho, gC = g, gR = g, grho = g, m,
               bC, aC, KC, bR, aR, KR, brho, arho1, Krho1, arho2, Krho2,
               regime = :S, target_mode = -1,
               label = "$label  (strength=$(round(strength; digits = 2)))")
end

"""
    build_edges(regime; fCrho, frhoC, fRC, frhoR, dC, dR, drho, ...)

Core builder: given the four target Jacobian edge values at u*=(1,1,1)
(all with regime-legal signs) and diffusion coefficients, return the
RDParams with matching Hill kinetics.  `regime==:P` uses an activating
Rac->Rho Hill term (frhoR must be > 0); `regime==:S` an inhibiting one
(frhoR must be < 0).
"""
function build_edges(regime::Symbol; fCrho, frhoC, fRC, frhoR,
        dC = 0.1, dR = 10.0, drho = 0.1, g = 1.0, m = 3,
        KC = 0.2, Krho1 = 0.2, KR = 5.0, Krho2 = 5.0, label = "")
    uC = uR = urho = 1.0
    dRterm = regime === :P ? dhillA(uR, Krho2, m) : dhillI(uR, Krho2, m)
    hRterm = regime === :P ? hillA(uR, Krho2, m)  : hillI(uR, Krho2, m)
    aC    = fCrho / dhillI(urho, KC,    m)
    arho1 = frhoC / dhillI(uC,   Krho1, m)
    aR    = fRC   / dhillA(uC,   KR,    m)
    arho2 = frhoR / dRterm
    bC   = g*uC   - aC   * hillI(urho, KC,    m)
    bR   = g*uR   - aR   * hillA(uC,   KR,    m)
    brho = g*urho - arho1* hillI(uC,   Krho1, m) - arho2 * hRterm
    RDParams(; dC, dR, drho, gC = g, gR = g, grho = g, m,
               bC, aC, KC, bR, aR, KR, brho, arho1, Krho1, arho2, Krho2,
               regime, target_mode = -1, label)
end

"""
    scan_regime(regime; ndraw, seed)

Random search over sign-legal parameter space of the given regime.
Returns counts: how many draws had a linearly stable homogeneous state,
and of those, how many admit a finite-wavelength (Turing) instability.
This is the numerical form of the Regime Dichotomy: ~0 for :S, > 0 for :P.
"""
function scan_regime(regime::Symbol; ndraw = 4000, seed = 2024, nmax = 20)
    rng = _rng(seed)
    nstable = 0; nturing = 0
    maxg_turing = Float64[]           # max finite-wavelength growth among stable-homog draws
    for _ in 1:ndraw
        lo, hi = 0.2, 6.0
        q = exp(log(lo) + rand(rng)*(log(hi)-log(lo)))
        r = exp(log(lo) + rand(rng)*(log(hi)-log(lo)))
        w = exp(log(lo) + rand(rng)*(log(hi)-log(lo)))
        dR = exp(log(2.0) + rand(rng)*(log(60.0)-log(2.0)))
        dslow = exp(log(0.01) + rand(rng)*(log(1.0)-log(0.01)))
        frhoR = regime === :P ? +w : -w
        p = build_edges(regime; fCrho = -q, frhoC = -q, fRC = r, frhoR = frhoR,
                        dC = dslow, dR = dR, drho = dslow)
        d = dispersion(p; nmax = nmax)
        if d.homog_stable
            nstable += 1
            g1 = maximum(d.growth[2:end])
            push!(maxg_turing, g1)
            d.turing && (nturing += 1)
        end
    end
    return (; regime, ndraw, nstable, nturing,
              frac_turing = nturing / max(nstable, 1), maxg_turing)
end

"Number of distinct positive homogeneous steady states (multi-start Newton) and how many are linearly stable (n=0)."
function count_homog_states(p::RDParams; grid = 0.05:0.35:4.0)
    found = Vector{Vector{Float64}}()
    for a in grid, b in grid, c in grid
        u, ok = homogeneous_state(p; u0 = [a, b, c])
        ok || continue
        all(u .> 0) || continue
        any(v -> norm(v .- u) < 1e-4, found) || push!(found, u)
    end
    nstable = 0
    for u in found
        ev = eigvals(reaction_jac(u, p))
        maximum(real, ev) < 1e-8 && (nstable += 1)
    end
    return (; nstates = length(found), nstable, states = found)
end

# =====================================================================
#  1-D PDE solver:  u_t = D u_xx + F(u)  on [0,2pi], no-flux
#  Method of lines + IMEX Euler (implicit diffusion, explicit reaction).
# =====================================================================
"""
    simulate(p; N, L, tmax, dt, seed, amp, tol)

Integrate to (numerical) steady state from a small random perturbation of
the homogeneous state.  Returns a NamedTuple with the grid `x`, the final
profiles `C,R,rho`, the homogeneous state `ustar`, the stopping time `t`,
and the final max |du/dt| `resid`.
"""
function simulate(p::RDParams; N = 257, L = 2π, tmax = 4000.0, dt = 4e-3,
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

# --------------------------- convenience ----------------------------
"Descriptive one-liner used in logging / CSV."
function summarise(p::RDParams, d)
    reg = p.regime === :P ? "PATTERN" : "SWITCH "
    @sprintf("%s | %-40s | homog_stable=%-5s | n*=%2d | growth*=%+.4f | turing=%s",
             reg, p.label, string(d.homog_stable), d.nstar, d.growth_star, string(d.turing))
end
