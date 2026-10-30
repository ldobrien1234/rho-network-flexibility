# =====================================================================
#  onset_sweeps.jl
#
#  Quasistatic sweeps across the first Turing bifurcation in cycle space.
#
#  For k = 1, 2, 3 the k..k region of Theorem A has an edge on the neutral
#  line of mode k,
#     T_k :  g_k(C2, C3) = prod_i (f_ii - d_i k^2) - C2 (f_RR - d_R k^2) + C3 = 0,
#  with mode k unstable where g_k > 0.  The path runs along the unit normal
#  of T_k through the midpoint of that edge,
#     (C2, C3)(s) = foot + s * nhat,      s < 0 stable,  s > 0 only mode k unstable,
#  so s is signed distance from the linear onset.  At every step the kinetics
#  are re-lifted (u* stays at (1,1,1)), the previous profile plus a tiny
#  perturbation is relaxed to steady state, and the mode-k amplitude is
#  recorded.  Forward (stable -> unstable) then backward.
#
#  Supercritical onset: amplitude grows from zero like sqrt(s), forward and backward agree.
#  Subcritical onset:   finite-amplitude jump at s >= 0, pattern persists to s < 0 on the way back.
#
#  The second half sweeps across the half-plane boundary W_0 (uniform Hopf).
#
#  Run:   julia -t 3 onset_sweeps.jl        Output in output/sweeps/.
# =====================================================================

include("cycle_atlas.jl")

const SWEEP_DIR = joinpath(@__DIR__, "output", "sweeps")

gk(C2, C3, k; f = F_DIAG, d = D_DIFF) = prod(f[i] - d[i]*k^2 for i in 1:3) - C2*(f[2] - d[2]*k^2) + C3

"""
Point at fraction `t` along the edge of the k..k region lying on T_k (ends ordered by C2; t = 0.5 is
the midpoint), and the unit normal pointing to g_k > 0.
"""
function onset_path(k; t = 0.5, R = theorem_regions())
    P = R[(k, k)]
    scale = maximum(abs, last.(P))
    on = filter(v -> abs(gk(v..., k)) < 1e-8scale, P)       # vertices on T_k (the edge may hold a collinear vertex)
    length(on) >= 2 || error("the $k..$k region has no edge on T_$k")
    a = on[argmin(first.(on))];  b = on[argmax(first.(on))]   # ends of the T_k edge
    foot = a .+ t .* (b .- a)
    gr = (-(F_DIAG[2] - D_DIFF[2]*k^2), 1.0)                 # gradient of g_k
    nhat = gr ./ hypot(gr...)
    return foot, nhat
end

"Unstable modes along the path; returns the s-range on which only mode k is unstable (s > 0)."
function path_extent(k, foot, nhat; smax = 5.0, ds = 1e-3)
    s = ds
    while s < smax
        q = foot .+ s .* nhat
        u, _ = unstable_modes(lift(q...))
        u == [k] || break
        s += ds
    end
    s_hi = s - ds
    s = -ds
    while s > -smax
        q = foot .+ s .* nhat
        u, _ = unstable_modes(lift(q...))
        isempty(u) || break
        s -= ds
    end
    return s + ds, s_hi                                      # (most negative stable s, largest single-mode s)
end

"LU factors of the implicit diffusion step (Neumann Laplacian on [0, L_DOM]), as in `simulate`."
function imex_factors(p::RDParams, N; dt = 0.01)
    h = L_DOM / (N - 1)
    dl = fill(1.0 / h^2, N - 1);  d0 = fill(-2.0 / h^2, N);  du = fill(1.0 / h^2, N - 1)
    du[1] = 2.0 / h^2;  dl[end] = 2.0 / h^2
    Lap = Tridiagonal(dl, d0, du);  I_N = Diagonal(ones(N))
    (lu(I_N - dt*p.dC*Lap), lu(I_N - dt*p.dR*Lap), lu(I_N - dt*p.drho*Lap))
end

evolve!(C, R, rho, p::RDParams, facs, dt, nsteps, tol) =
    _evolve!(C, R, rho, facs, true, p.bC, p.aC, p.KC, p.gC, p.bR, p.aR, p.KR, p.gR,
             p.brho, p.arho1, p.Krho1, p.arho2, p.Krho2, p.grho, p.m, dt, nsteps, tol)

"Relax (C,R,rho) in place under parameters p until max|du/dt| < tol or t > tmax."
relax!(C, R, rho, p::RDParams; dt = 0.01, tmax = 5000.0, tol = 1e-9) =
    evolve!(C, R, rho, p, imex_factors(p, length(C); dt), dt, ceil(Int, tmax/dt), tol)

"""
    sweep(k; nstep, N, kick, seed)

Forward then backward quasistatic sweep across T_k.  Returns a vector of
NamedTuples (dir, s, C2, C3, growth_k, a_k, nd, share, resid, t).
"""
function sweep(k; t = 0.5, nstep = 70, N = 257, kick = 1e-6, seed = 1, s_min = -1.0, frac_hi = 0.85)
    foot, nhat = onset_path(k; t)
    s_lo, s_hi = path_extent(k, foot, nhat)
    ss = collect(range(max(s_min, s_lo), frac_hi * s_hi; length = nstep))
    x = range(0, L_DOM; length = N)
    rng = Xoshiro(seed)
    C = ones(N);  R = ones(N);  rho = ones(N)
    out = []
    for (dir, path) in ((:forward, ss), (:backward, reverse(ss)))
        for s in path
            q = foot .+ s .* nhat
            p = lift(q...)
            C .+= kick .* randn(rng, N);  R .+= kick .* randn(rng, N);  rho .+= kick .* randn(rng, N)
            # A kick of 1e-6 growing at rate lam keeps |du/dt| below tol for a long time,
            # so run at least 12/|lam| (growth by e^12) before the tolerance may stop it.
            d = dispersion(p; nmax = 12, ustar = [1.0, 1.0, 1.0])
            lam = d.growth[k + 1]
            tmin = min(12 / max(abs(lam), 1e-4), 40000.0)
            t1, _ = relax!(C, R, rho, p; tmax = tmin, tol = 0.0)
            t2, resid = relax!(C, R, rho, p; tmax = 40000.0)
            t = t1 + t2
            a = cosine_amps(C, x; nmax = 12)
            amp = maximum(maximum(u) - minimum(u) for u in (C, R, rho))
            state = amp > 1e-4 ? :pattern : abs(sum(C)/N - 1) < 1e-3 ? :ustar : :other_uniform
            nd = state === :pattern ? argmax(abs.(a)) : 0
            push!(out, (; k, dir, s, C2 = q[1], C3 = q[2], growth_k = lam, a_k = a[k], nd, state,
                          meanC = sum(C)/N, meanrho = sum(rho)/N,
                          share = nd == 0 ? NaN : a[nd]^2 / sum(abs2, a), resid, t,
                          C = copy(C), R = copy(R), rho = copy(rho)))
        end
    end
    return out, (; foot, nhat, s_lo, s_hi)
end

# =====================================================================
#  Output
# =====================================================================
function write_sweep_csv(all)
    mkpath(SWEEP_DIR)
    open(joinpath(SWEEP_DIR, "sweeps.csv"), "w") do io
        println(io, "k,dir,s,C2,C3,growth_k,a_k,state,dominant_mode,share,meanC,meanrho,resid,t")
        for r in all
            @printf(io, "%d,%s,%.6g,%.8g,%.8g,%.5g,%.5g,%s,%d,%.4g,%.5g,%.5g,%.3g,%.1f\n",
                    r.k, r.dir, r.s, r.C2, r.C3, r.growth_k, r.a_k, r.state, r.nd, r.share,
                    r.meanC, r.meanrho, r.resid, r.t)
        end
    end
end

stateclass(r) = r.state === :ustar ? -1 : r.nd                # -1 u*, 0 other uniform, n pattern mode
statelabel(r) = r.state === :ustar ? "u*" : r.state === :other_uniform ? "other uniform" : "mode $(r.nd)"

ampC(r) = maximum(r.C) - minimum(r.C)

"Pattern amplitude (max - min of Cdc42) along each sweep, forward filled / backward open, coloured by state."
function plot_sweeps(results)
    mkpath(SWEEP_DIR)
    panels = []
    for (rows, info) in results
        k = rows[1].k
        fw = filter(r -> r.dir === :forward, rows);  bw = filter(r -> r.dir === :backward, rows)
        pa = plot(; xlabel = "s  (signed distance past T_$k in (C₂,C₃))", ylabel = "max C − min C",
                  title = "sweep across T_$k  (only mode $k linearly unstable for s > 0):  forward ●  backward ○",
                  titlefontsize = 9, legend = :topleft)
        vline!(pa, [0.0]; color = :black, ls = :dash, label = "linear onset (T_$k)")
        plot!(pa, [r.s for r in fw], ampC.(fw); color = :gray50, lw = 1, label = "")
        plot!(pa, [r.s for r in bw], ampC.(bw); color = :gray50, lw = 1, ls = :dot, label = "")
        for cls in sort(unique(stateclass(r) for r in rows))
            f_ = filter(r -> stateclass(r) == cls, fw);  b_ = filter(r -> stateclass(r) == cls, bw)
            col = cls == -1 ? :black : modecol(cls)
            lab = cls == -1 ? "u* (uniform)" : cls == 0 ? "other uniform state (Rho-high)" : "mode $cls pattern"
            scatter!(pa, [r.s for r in f_], ampC.(f_); color = col, ms = 4, msw = 0, label = lab)
            scatter!(pa, [r.s for r in b_], ampC.(b_); color = :white, msc = col, ms = 5, msw = 1.2, label = "")
        end
        push!(panels, pa)
    end
    plt = plot(panels...; layout = (length(results), 1), size = (900, 300 * length(results)), left_margin = 4Plots.mm)
    savefig(plt, joinpath(SWEEP_DIR, "onset_sweeps.png"))
    plt
end

"Profiles per sweep: first state past onset (forward), end of forward sweep, last pattern on the way back."
function plot_sweep_profiles(results)
    panels = []
    for (rows, info) in results
        fw = filter(r -> r.dir === :forward, rows);  bw = filter(r -> r.dir === :backward, rows)
        i_on = something(findfirst(r -> r.state !== :ustar, fw), length(fw))
        i_bk = something(findlast(r -> r.state === :pattern, bw), 1)
        for (r, what) in ((fw[i_on], "forward, first past onset"), (fw[end], "forward, end"), (bw[i_bk], "backward, last pattern"))
            x = range(0, L_DOM; length = length(r.C))
            plt = plot(x, r.C; color = CCOL, lw = 2, label = "Cdc42", legend = isempty(panels) ? :topright : false,
                       title = @sprintf("k = %d, %s: s = %.3g → %s", r.k, what, r.s,
                                        r.state === :pattern ? "mode $(r.nd)" : r.state === :ustar ? "u*" : "other uniform"),
                       titlefontsize = 8, xlims = (0, π), xticks = ([0, π/2, π], ["0", "π/2", "π"]))
            plot!(plt, x, r.R; color = RCOL, lw = 2, label = "Rac")
            plot!(plt, x, r.rho; color = RHOCOL, lw = 2, label = "Rho")
            push!(panels, plt)
        end
    end
    plt = plot(panels...; layout = (length(results), 3), size = (1200, 260 * length(results)))
    savefig(plt, joinpath(SWEEP_DIR, "onset_profiles.png"))
    plt
end

# =====================================================================
#  Hopf crossing of the half-plane boundary
#
#  W_0 :  C3 = a C2 + b,  a = -(f_CC + f_rhorho),  b = (f_CC+f_RR)(f_CC+f_rhorho)(f_RR+f_rhorho)
#  is the Routh-Hurwitz condition a1 a2 = a3 for J itself.  On it J has
#  eigenvalues +-i omega, omega = sqrt(a2); outside the half-plane the uniform
#  mode k = 0 is an unstable complex pair, so uniform oscillations should
#  start with period near 2 pi / omega.  For C2 below ~34 (canonical set)
#  every k >= 1 mode stays stable on both sides, so the crossing is a pure
#  k = 0 Hopf.  Oscillations never reach a steady state, so each step
#  integrates a settling time and then records the spatial means over ~12
#  periods.
# =====================================================================
# The default IMEX Euler step treats the reaction explicitly: for eigenvalues +-i omega it
# amplifies by |1 + i omega dt| > 1, an artificial growth rate ~ omega^2 dt / 2 (~0.01 here)
# that moves the Hopf point.  Steady states and real eigenvalues are unaffected, but the
# Hopf sweep needs the reaction stepped accurately, so it uses RK4 for the reaction
# followed by an implicit diffusion step (Lie splitting; exact to RK4 order for uniform states).
function react!(dC, dR, drho, C, R, rho, p::RDParams)
    @inbounds for i in eachindex(C)
        dC[i]   = p.bC + p.aC * hillI(rho[i], p.KC, p.m) - p.gC * C[i]
        dR[i]   = p.bR + p.aR * hillA(C[i], p.KR, p.m) - p.gR * R[i]
        drho[i] = p.brho + p.arho1 * hillI(C[i], p.Krho1, p.m) + p.arho2 * hillA(R[i], p.Krho2, p.m) - p.grho * rho[i]
    end
end

"`nsteps` steps of RK4 reaction + implicit diffusion, in place."
function evolve_rk4!(C, R, rho, p::RDParams, facs, dt, nsteps)
    N = length(C)
    k = [zeros(N) for _ in 1:12];  tC = zeros(N);  tR = zeros(N);  trho = zeros(N)
    for _ in 1:nsteps
        react!(k[1], k[2], k[3], C, R, rho, p)
        @. tC = C + dt/2*k[1];  @. tR = R + dt/2*k[2];  @. trho = rho + dt/2*k[3]
        react!(k[4], k[5], k[6], tC, tR, trho, p)
        @. tC = C + dt/2*k[4];  @. tR = R + dt/2*k[5];  @. trho = rho + dt/2*k[6]
        react!(k[7], k[8], k[9], tC, tR, trho, p)
        @. tC = C + dt*k[7];  @. tR = R + dt*k[8];  @. trho = rho + dt*k[9]
        react!(k[10], k[11], k[12], tC, tR, trho, p)
        @. C   += dt/6 * (k[1] + 2k[4] + 2k[7] + k[10])
        @. R   += dt/6 * (k[2] + 2k[5] + 2k[8] + k[11])
        @. rho += dt/6 * (k[3] + 2k[6] + 2k[9] + k[12])
        ldiv!(facs[1], C);  ldiv!(facs[2], R);  ldiv!(facs[3], rho)
    end
end

"Foot on W_0 at the given C2, and the unit normal pointing out of the half-plane."
function hopf_path(C2; f = F_DIAG)
    a = -(f[1] + f[3]);  b = (f[1] + f[2]) * (f[1] + f[3]) * (f[2] + f[3])
    return (C2, a*C2 + b), (a, -1.0) ./ hypot(a, 1.0)
end

"Largest s > 0 on which only k = 0 is unstable, as a complex pair; errors if s < 0 is not stable down to s_min."
function hopf_extent(foot, nhat; s_min = -0.5, smax = 3.0, ds = 1e-3)
    for s in range(s_min, -ds; step = ds)
        u, _ = unstable_modes(lift((foot .+ s .* nhat)...))
        isempty(u) || error("mode(s) $u unstable at s = $s on the stable side of W_0")
    end
    s = ds
    while s < smax
        u, d = unstable_modes(lift((foot .+ s .* nhat)...))
        (u == [0] && d.complex[1]) || break
        s += ds
    end
    return s - ds
end

"Integrate for time T, sampling the spatial means of C, R, rho and the spatial range of C every dts."
function record!(C, R, rho, p::RDParams; T, dts, dt = 0.01)
    N = length(C);  facs = imex_factors(p, N; dt)
    nsub = max(1, round(Int, dts / dt));  ns = round(Int, T / (nsub * dt))
    tt = (1:ns) .* (nsub * dt)
    mC = zeros(ns);  mR = zeros(ns);  mrho = zeros(ns);  spread = zeros(ns)
    for i in 1:ns
        evolve_rk4!(C, R, rho, p, facs, dt, nsub)
        mC[i] = sum(C) / N;  mR[i] = sum(R) / N;  mrho[i] = sum(rho) / N
        spread[i] = maximum(C) - minimum(C)
    end
    return (; tt = collect(tt), mC, mR, mrho, spread)
end

"Mean period from upward zero crossings of v - mean(v) (linear interpolation); NaN if fewer than 3."
function period_of(tt, v)
    w = v .- sum(v) / length(v)
    ups = [tt[i] - w[i] * (tt[i+1] - tt[i]) / (w[i+1] - w[i]) for i in 1:length(w)-1 if w[i] < 0 <= w[i+1]]
    length(ups) < 3 ? NaN : (ups[end] - ups[1]) / (length(ups) - 1)
end

"""
    hopf_sweep(; C2, nstep, s_min, frac_hi, N, kick, seed)

Forward then backward sweep across W_0 at the given C2.  Each step: settle for
max(40 periods, 15/|Re lambda_0|) (capped), then record ~12 periods.
"""
function hopf_sweep(; C2 = 33.0, nstep = 40, s_min = -0.5, frac_hi = 0.85, N = 65, kick = 1e-6, seed = 1)
    foot, nhat = hopf_path(C2)
    s_hi = hopf_extent(foot, nhat; s_min)
    ss = collect(range(s_min, frac_hi * s_hi; length = nstep))
    rng = Xoshiro(seed)
    C = ones(N);  R = ones(N);  rho = ones(N)
    out = []
    for (dir, path) in ((:forward, ss), (:backward, reverse(ss)))
        for s in path
            q = foot .+ s .* nhat
            p = lift(q...)
            ev = eigvals(reaction_jac([1.0, 1.0, 1.0], p))
            i0 = argmax(real.(ev));  lam = real(ev[i0]);  omega = abs(imag(ev[i0]))
            Tlin = 2π / omega
            C .+= kick .* randn(rng, N);  R .+= kick .* randn(rng, N);  rho .+= kick .* randn(rng, N)
            Tset = clamp(15 / max(abs(lam), 1e-4), 40Tlin, 8000.0)
            evolve_rk4!(C, R, rho, p, imex_factors(p, N; dt = 0.01), 0.01, round(Int, Tset / 0.01))
            rec = record!(C, R, rho, p; T = 12Tlin, dts = Tlin / 50)
            amp = (maximum(rec.mC) - minimum(rec.mC)) / 2
            # a limit cycle keeps its amplitude across the record; a slowly decaying transient does not
            h = length(rec.mC) ÷ 2
            amp1 = maximum(rec.mC[1:h]) - minimum(rec.mC[1:h]);  amp2 = maximum(rec.mC[h+1:end]) - minimum(rec.mC[h+1:end])
            state = amp > 1e-6 && amp2 > 0.98 * amp1 ? :oscillation :
                    amp > 1e-6 && abs(rec.mC[end] - 1) < 1e-2 ? :ustar :
                    abs(rec.mC[end] - 1) < 1e-3 ? :ustar : :other_uniform
            push!(out, (; dir, s, C2 = q[1], C3 = q[2], lam, omega, Tlin, amp,
                          period = state === :oscillation ? period_of(rec.tt, rec.mC) : NaN,
                          maxspread = maximum(rec.spread), state, rec))
        end
    end
    return out, (; foot, nhat, s_hi)
end

function write_hopf_csv(rows)
    mkpath(SWEEP_DIR)
    open(joinpath(SWEEP_DIR, "hopf_sweep.csv"), "w") do io
        println(io, "dir,s,C2,C3,re_lambda0,omega0,period_linear,state,amp_meanC,period,max_spatial_range_C")
        for r in rows
            @printf(io, "%s,%.6g,%.8g,%.8g,%.5g,%.5g,%.5g,%s,%.5g,%.5g,%.3g\n", r.dir, r.s, r.C2, r.C3,
                    r.lam, r.omega, r.Tlin, r.state, r.amp, r.period, r.maxspread)
        end
    end
end

"Least-squares amp^2 = c (s - s0) over the first `n` oscillating forward steps."
function fit_hopf(rows; n = 6)
    fw = filter(r -> r.dir === :forward && r.state === :oscillation, rows)
    length(fw) < 3 && return nothing
    fw = fw[1:min(n, end)]
    b = [ones(length(fw)) [r.s for r in fw]] \ [r.amp^2 for r in fw]
    (; c = b[2], s0 = -b[1] / b[2])
end

function plot_hopf(rows, info)
    fw = filter(r -> r.dir === :forward, rows);  bw = filter(r -> r.dir === :backward, rows)
    osc(v) = filter(r -> r.state === :oscillation, v)
    pa = plot(; xlabel = "s  (signed distance past W₀ in (C₂,C₃))", ylabel = "oscillation amplitude of mean C",
              title = @sprintf("Hopf sweep across the half-plane boundary at C₂ ≈ %.3g:  forward ●  backward ○", info.foot[1]),
              titlefontsize = 9, legend = :topleft)
    vline!(pa, [0.0]; color = :black, ls = :dash, label = "linear onset (W₀)")
    scatter!(pa, [r.s for r in fw], [r.amp for r in fw]; color = :black, ms = 4, msw = 0, label = "forward")
    scatter!(pa, [r.s for r in bw], [r.amp for r in bw]; color = :white, msc = :black, ms = 5, label = "backward")
    fit = fit_hopf(rows)
    if fit !== nothing
        xs = range(max(fit.s0, 0.0), maximum(r.s for r in rows); length = 100)
        plot!(pa, xs, sqrt.(max.(fit.c .* (xs .- fit.s0), 0.0)); color = :red,
              label = @sprintf("√(%.3g (s − %.2g))", fit.c, fit.s0))
    end
    pb = plot(; xlabel = "s", ylabel = "period", title = "period: measured vs 2π/ω of the linearisation",
              titlefontsize = 9, legend = :topright)
    plot!(pb, [r.s for r in fw], [r.Tlin for r in fw]; color = :gray50, lw = 2, label = "2π/ω (linear, at u*)")
    scatter!(pb, [r.s for r in osc(fw)], [r.period for r in osc(fw)]; color = :black, ms = 4, msw = 0, label = "measured, forward")
    scatter!(pb, [r.s for r in osc(bw)], [r.period for r in osc(bw)]; color = :white, msc = :black, ms = 5, label = "measured, backward")
    vline!(pb, [0.0]; color = :black, ls = :dash, label = "")
    ofw = osc(fw)
    picks = isempty(ofw) ? [] : ofw[unique(round.(Int, range(1, length(ofw); length = 3)))]
    pcs = map(picks) do r
        t0 = r.rec.tt .- r.rec.tt[1]
        plt = plot(t0, r.rec.mC; color = CCOL, lw = 2, label = "⟨Cdc42⟩", titlefontsize = 8,
                   title = @sprintf("s = %.3g:  max spatial range of C = %.1e", r.s, r.maxspread),
                   xlabel = "t", legend = r === picks[1] ? :topright : false)
        plot!(plt, t0, r.rec.mR; color = RCOL, lw = 2, label = "⟨Rac⟩")
        plot!(plt, t0, r.rec.mrho; color = RHOCOL, lw = 2, label = "⟨Rho⟩")
    end
    top = plot(pa, pb; layout = (1, 2))
    plt = isempty(pcs) ? plot(top; size = (1200, 380)) :
        plot(top, plot(pcs...; layout = (1, length(pcs))); layout = grid(2, 1; heights = [0.55, 0.45]), size = (1200, 700),
             left_margin = 4Plots.mm, bottom_margin = 3Plots.mm)
    savefig(plt, joinpath(SWEEP_DIR, "hopf_sweep.png"))
    plt
end

function report_hopf(rows, info)
    fw = filter(r -> r.dir === :forward, rows);  bw = filter(r -> r.dir === :backward, rows)
    @printf("Hopf: foot (C2,C3) = (%.4f, %.4f), normal = (%.3f, %.3f), only k = 0 (complex) unstable up to s = %.3g\n",
            info.foot..., info.nhat..., info.s_hi)
    on  = findfirst(r -> r.state === :oscillation, fw)
    off = findfirst(r -> r.state !== :oscillation, bw)
    @printf("   forward: oscillation from s = %s;  backward: oscillation lost at s = %s  (step %.3g)\n",
            on === nothing ? "-" : @sprintf("%.4g", fw[on].s), off === nothing ? "never" : @sprintf("%.4g", bw[off].s),
            fw[2].s - fw[1].s)
    @printf("   states forward: %s;  backward: %s\n", join(unique(string(r.state) for r in fw), " → "),
            join(unique(string(r.state) for r in bw), " → "))
    fit = fit_hopf(rows)
    fit === nothing || @printf("   near-onset fit amp² = %.4g (s − %.3g)\n", fit.c, fit.s0)
    for r in filter(r -> r.state === :oscillation, fw)[1:min(3, end)]
        @printf("   s = %.3g: amp %.3g, period %.4g (linear 2π/ω = %.4g), max spatial range of C %.1e\n",
                r.s, r.amp, r.period, r.Tlin, r.maxspread)
    end
end

# =====================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    ks = [1, 2, 3]
    results = Vector{Any}(undef, length(ks))
    @threads for i in eachindex(ks)
        results[i] = sweep(ks[i])
    end
    for (rows, info) in results
        k = rows[1].k
        @printf("k = %d: foot (C2,C3) = (%.4f, %.4f), normal = (%.3f, %.3f), stable down to s = %.3g, only mode %d up to s = %.3g\n",
                k, info.foot..., info.nhat..., info.s_lo, k, info.s_hi)
        fw = filter(r -> r.dir === :forward, rows);  bw = filter(r -> r.dir === :backward, rows)
        on  = findfirst(r -> r.state !== :ustar, fw)
        off = findfirst(r -> r.state === :ustar, bw)
        ds = fw[2].s - fw[1].s
        @printf("   forward: leaves u* at s = %s;  backward: back on u* at s = %s  (step %.3g)\n",
                on === nothing ? "-" : @sprintf("%.4g", fw[on].s),
                off === nothing ? "never (down to s_min)" : @sprintf("%.4g", bw[off].s), ds)
        @printf("   forward states: %s\n", join(unique(statelabel(r) for r in fw), " → "))
        @printf("   backward states: %s\n", join(unique(statelabel(r) for r in bw), " → "))
        i_on = findfirst(r -> r.state !== :ustar, fw)
        i_on === nothing || @printf("   amplitude (max-min C) at first step past onset: %.3g;  forward end: %.3g\n",
                                    ampC(fw[i_on]), ampC(fw[end]))
        nslow = count(r -> r.resid > 1e-7, rows)
        nslow > 0 && @printf("   %d steps not fully relaxed (resid > 1e-7)\n", nslow)
    end
    write_sweep_csv(vcat(first.(results)...))
    plot_sweeps(results)
    plot_sweep_profiles(results)

    hrows, hinfo = hopf_sweep()
    report_hopf(hrows, hinfo)
    write_hopf_csv(hrows)
    plot_hopf(hrows, hinfo)
    println("figures and CSV in ", SWEEP_DIR)
end
