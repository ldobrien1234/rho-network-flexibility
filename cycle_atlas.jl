# =====================================================================
#  cycle_atlas.jl
#
#  Attractor atlas over cycle space (C2, C3).
#
#  C2 = f_Crho f_rhoC   (Cdc42 <-> Rho 2-cycle)
#  C3 = f_RC f_rhoR f_Crho   (Cdc42 -> Rac -> Rho -| Cdc42 3-cycle)
#
#  Theorem A (plot_parameter_region.m): on polygon ∩ half-plane the unstable
#  modes of J - k^2 D are exactly k = kappa1..kappa2 (needs tr J < 0 too).
#  For each window we sample points in that set, lift each point to Hill
#  kinetics with the prescribed Jacobian at u* = (1,1,1), check the linear
#  window, and integrate the PDE from several initial conditions to see
#  which patterns the nonlinear system actually settles into.
#
#  Domain: x in [0, pi] with no-flux boundaries, so the Laplacian modes are
#  exactly cos(n x), B = -n^2, matching the theorem.
#
#  Run:   julia -t auto cycle_atlas.jl
#  Output in output/atlas/.
# =====================================================================

include("cdc42_rac_rho.jl")
using Plots
using Base.Threads

# --------------------------- canonical parameter set -----------------
const F_DIAG = (-5.0, -1.0, -5.0)       # (f_CC, f_RR, f_rhorho)
const D_DIFF = (0.05, 0.5, 0.05)        # (d_C, d_R, d_rho)
const L_DOM  = π                        # Neumann on [0, pi]: modes cos(n x)
const KMAX   = 4                        # windows kappa1..kappa2 with kappa2 <= KMAX
const OUTDIR = joinpath(@__DIR__, "output", "atlas")
const CCOL = "#1b7837"; const RCOL = "#2166ac"; const RHOCOL = "#b2182b"   # as in run_simulations.jl

# =====================================================================
#  Theorem A regions (port of plot_parameter_region.m)
# =====================================================================
"Vertex (x,y) pairs of the multiple-unstable-mode table."
function vertex_pairs(S, m, sn, s1, M)
    S <= m + sn  && return Tuple{Float64,Float64}[]
    S <  m + M   && return [(m, sn), (m, 0.0), (sn, 0.0)]
    S <  s1 + sn && return [(m, M), (m, sn), (M, 0.0), (s1, 0.0)]
    S <  s1 + M  && return [(m, M), (m, sn), (M, 0.0), (sn, 0.0)]
    return [(s1, M), (s1, sn), (m, M), (m, sn)]
end

function sort_ccw(P)
    c = (sum(first, P) / length(P), sum(last, P) / length(P))
    sort(P; by = v -> atan(v[2] - c[2], v[1] - c[1]))
end

"Clip a convex polygon to the half-plane a*x + b*y + c >= 0 (Sutherland-Hodgman)."
function clip(P, a, b, c)
    out = Tuple{Float64,Float64}[]
    n = length(P)
    for i in 1:n
        p = P[i]; q = P[mod1(i + 1, n)]
        sp = a*p[1] + b*p[2] + c;  sq = a*q[1] + b*q[2] + c
        sp >= 0 && push!(out, p)
        if sp * sq < 0
            t = sp / (sp - sq)
            push!(out, (p[1] + t*(q[1] - p[1]), p[2] + t*(q[2] - p[2])))
        end
    end
    return out
end

polyarea(P) = length(P) < 3 ? 0.0 :
    0.5 * abs(sum(P[i][1]*P[mod1(i+1, length(P))][2] - P[mod1(i+1, length(P))][1]*P[i][2] for i in eachindex(P)))

function inpoly(P, q)
    n = length(P)
    all((P[mod1(i+1, n)][1] - P[i][1])*(q[2] - P[i][2]) -
        (P[mod1(i+1, n)][2] - P[i][2])*(q[1] - P[i][1]) >= 0 for i in 1:n)
end

"""
    theorem_regions(f, d; kmax)

Theorem A region for every window (a,b), 1 <= a <= b <= kmax, intersected
with the half-plane and with {C2 > 0, C3 < 0} (regime P with mutual
C <-> rho inhibition).  Returns Dict((a,b) => ccw polygon).
"""
function theorem_regions(f = F_DIAG, d = D_DIFF; kmax = KMAX)
    phi = collect(f) ./ collect(d)
    S   = -sum(phi)
    xi2 = phi[1]*phi[2] + phi[1]*phi[3] + phi[2]*phi[3]
    xi3 = prod(phi)
    vert(x, y) = (z = S - x - y;
                  C2 = d[1]*d[3] * (xi2 - x*y - (x + y)*z);
                  C3 = prod(d) * (phi[2]*C2/(d[1]*d[3]) - x*y*z - xi3);
                  (C2, C3))
    ah = -(f[1] + f[3])
    bh = (f[1] + f[2]) * (f[1] + f[3]) * (f[2] + f[3])
    R = Dict{Tuple{Int,Int},Vector{Tuple{Float64,Float64}}}()
    for a in 1:kmax, b in a:kmax
        m  = -(b + 1)^2;  s1 = -b^2;  sn = -a^2;  M = a == 1 ? 0.0 : -(a - 1)^2
        pr = vertex_pairs(S, m, sn, s1, M)
        isempty(pr) && continue
        P = sort_ccw([vert(p...) for p in pr])
        P = clip(P, -ah, 1.0, -bh)       # half-plane  C3 > ah*C2 + bh
        P = clip(P, 1.0, 0.0, 0.0)       # C2 >= 0
        P = clip(P, 0.0, -1.0, 0.0)      # C3 <= 0
        polyarea(P) > 0 && (R[(a, b)] = P)
    end
    return R
end

"Centroid plus uniform random interior points of polygon P."
function sample_polygon(P, n; rng)
    c = (sum(first, P) / length(P), sum(last, P) / length(P))
    pts = [c]
    lo = (minimum(first, P), minimum(last, P));  hi = (maximum(first, P), maximum(last, P))
    while length(pts) < n
        q = (lo[1] + rand(rng)*(hi[1] - lo[1]), lo[2] + rand(rng)*(hi[2] - lo[2]))
        inpoly(P, q) && push!(pts, q)
    end
    return pts
end

# =====================================================================
#  Lift (C2, C3) to Hill kinetics with that Jacobian at u* = (1,1,1)
# =====================================================================
"""
    lift(C2, C3; fRC, m, KC, Krho1, KR, Krho2)

Edges: f_Crho = f_rhoC = -sqrt(C2), f_RC fixed, f_rhoR = C3/(f_RC f_Crho) > 0.
Diagonal f_ii = -g_i from F_DIAG (linear turnover).  Every basal rate b_i
must be >= 0 so production stays non-negative; errors otherwise.
"""
function lift(C2, C3; f = F_DIAG, d = D_DIFF, fRC = 1.0, m = 3,
              KC = 0.5, Krho1 = 0.5, KR = 2.0, Krho2 = 2.0, label = "")
    C2 > 0 && C3 < 0 || error("regime P needs C2 > 0 and C3 < 0 (got $C2, $C3)")
    fCrho = -sqrt(C2);  frhoC = -sqrt(C2)
    frhoR = C3 / (fRC * fCrho)
    gC, gR, grho = -f[1], -f[2], -f[3]
    aC    = fCrho / dhillI(1.0, KC,    m)
    arho1 = frhoC / dhillI(1.0, Krho1, m)
    aR    = fRC   / dhillA(1.0, KR,    m)
    arho2 = frhoR / dhillA(1.0, Krho2, m)
    bC   = gC   - aC * hillI(1.0, KC, m)
    bR   = gR   - aR * hillA(1.0, KR, m)
    brho = grho - arho1 * hillI(1.0, Krho1, m) - arho2 * hillA(1.0, Krho2, m)
    min(bC, bR, brho) >= 0 ||
        error(@sprintf("negative basal rate at (C2,C3)=(%.3g,%.3g): bC=%.3g bR=%.3g brho=%.3g", C2, C3, bC, bR, brho))
    RDParams(; dC = d[1], dR = d[2], drho = d[3], gC, gR, grho, m,
               bC, aC, KC, bR, aR, KR, brho, arho1, Krho1, arho2, Krho2,
               regime = :P, label)
end

"Unstable admissible modes 0..nmax of the lifted model at u* = (1,1,1)."
function unstable_modes(p::RDParams; nmax = 30)
    d = dispersion(p; nmax, ustar = [1.0, 1.0, 1.0])
    [n for (n, g) in zip(d.n, d.growth) if g > 1e-10], d
end

# =====================================================================
#  Classification of a final state on [0, pi]
# =====================================================================
"Cosine amplitudes a_n of v on [0, L] (trapezoid rule), n = 1..nmax."
function cosine_amps(v, x; nmax = 16)
    L = x[end] - x[1];  w = fill(step(x), length(x));  w[1] /= 2;  w[end] /= 2
    vm = v .- sum(w .* v) / L
    [2 / L * sum(w .* vm .* cos.(n .* x)) for n in 1:nmax]
end

"Local maxima of v (boundary points count, no-flux) rising at least `prom` of the range above the minimum."
function count_peaks(v; prom = 0.25)
    lo = minimum(v);  thr = lo + prom * (maximum(v) - lo)
    n = length(v);  c = 0
    v[1] >= v[2] && v[1] > thr && (c += 1)
    v[n] >= v[n-1] && v[n] > thr && (c += 1)
    for i in 2:n-1
        v[i] > v[i-1] && v[i] >= v[i+1] && v[i] > thr && (c += 1)
    end
    return c
end

"""
Summary of one simulation: kind (:uniform_ustar, :uniform_other, :pattern,
:unconverged), dominant cosine mode of C, its share of the cosine power,
peak count, amplitude, and C-rho / C-R correlations.
"""
function classify(r; x, amp_tol = 1e-4, tol = 1e-7)
    amp = maximum(maximum(u) - minimum(u) for u in (r.C, r.R, r.rho))
    corr(a, b) = (a = a .- sum(a)/length(a); b = b .- sum(b)/length(b);
                  s = sqrt(sum(abs2, a) * sum(abs2, b)); s > 0 ? sum(a .* b) / s : NaN)
    if r.resid > 100tol
        kind = :unconverged
    elseif amp < amp_tol
        kind = maximum(abs, [r.C[1], r.R[1], r.rho[1]] .- 1.0) < 1e-3 ? :uniform_ustar : :uniform_other
    else
        kind = :pattern
    end
    a = cosine_amps(r.C, x)
    nd = kind === :pattern ? argmax(abs.(a)) : 0
    share = kind === :pattern ? a[nd]^2 / sum(abs2, a) : NaN
    (; kind, nd, share, peaks = kind === :pattern ? count_peaks(r.C) : 0, amp,
       corr_Crho = corr(r.C, r.rho), corr_CR = corr(r.C, r.R), resid = r.resid, t = r.t)
end

"Distinct final states at one point, identified up to the reflection x -> L - x."
function distinct_states(finals; rtol = 1e-2)
    reps = Vector{Float64}[]
    for v in finals
        scale = max(maximum(abs, v), 1e-12)
        any(min(maximum(abs, v .- u), maximum(abs, reverse(v) .- u)) < rtol * scale for u in reps) ||
            push!(reps, v)
    end
    return length(reps)
end

# =====================================================================
#  Atlas driver
# =====================================================================
"""
    run_atlas(; npts, nnoise, N, dt, tmax, seed)

For every window: `npts` points (centroid first), linear check, then
`nnoise` random-noise ICs and one cosine-seeded IC per mode in the window.
Returns (runs, points) as vectors of NamedTuples.
"""
function run_atlas(; npts = 30, nnoise = 3, N = 257, dt = 0.01, tmax = 10000.0,
                   seed = 1, amp_noise = 1e-2, amp_mode = 5e-2, verbose = true)
    R = theorem_regions()
    rng = Xoshiro(seed)
    pts = []                                    # (window, C2, C3)
    for w in sort(collect(keys(R)))
        for q in sample_polygon(R[w], npts; rng)
            push!(pts, (w, q[1], q[2]))
        end
    end
    # jobs: (point index, ic symbol, seed / mode)
    jobs = Tuple{Int,Symbol,Int}[]
    for (i, (w, _, _)) in enumerate(pts)
        for s in 1:nnoise; push!(jobs, (i, :perturb, s)); end
        for k in w[1]:w[2]; push!(jobs, (i, :mode, k)); end
    end
    verbose && @printf("%d windows, %d points, %d simulations on %d threads\n",
                       length(R), length(pts), length(jobs), nthreads())

    params = [lift(C2, C3; label = @sprintf("w=%d..%d C2=%.4g C3=%.4g", w..., C2, C3)) for (w, C2, C3) in pts]
    lin = [unstable_modes(p) for p in params]
    x = range(0, L_DOM; length = N)

    results = Vector{Any}(undef, length(jobs))
    finals  = Vector{Vector{Float64}}(undef, length(jobs))
    done = Atomic{Int}(0)
    @threads :dynamic for j in eachindex(jobs)
        i, ic, s = jobs[j]
        r = ic === :perturb ?
            simulate(params[i]; N, L = L_DOM, dt, tmax, seed = s, amp = amp_noise, ic = :perturb) :
            simulate(params[i]; N, L = L_DOM, dt, tmax, seed = 1, amp = amp_mode, ic = :mode, seed_mode = s)
        results[j] = classify(r; x)
        finals[j]  = vcat(r.C, r.R, r.rho)
        nd = atomic_add!(done, 1) + 1
        verbose && nd % 50 == 0 && @printf("  %d / %d\n", nd, length(jobs))
    end

    runs = [(; point = i, window = pts[i][1], C2 = pts[i][2], C3 = pts[i][3], ic, icarg = s, results[j]...)
            for (j, (i, ic, s)) in enumerate(jobs)]
    points = map(eachindex(pts)) do i
        js = findall(jb -> jb[1] == i, jobs)
        unst, d = lin[i]
        w = pts[i][1]
        noise_modes = [results[j].nd for j in js if jobs[j][2] === :perturb]
        (; point = i, window = w, C2 = pts[i][2], C3 = pts[i][3],
           linear_ok = unst == collect(w[1]:w[2]),
           nstar = d.nstar, growth_star = d.growth_star,
           noise_modes, all_modes = sort(unique(results[j].nd for j in js)),
           nclasses = length(unique(results[j].nd for j in js)),   # mode n, or 0 = other uniform state
           nshapes = distinct_states(finals[js]),                    # also separates half-wavelength shifts
           any_unconverged = any(results[j].kind === :unconverged for j in js))
    end
    return runs, points, R, finals, x
end

# =====================================================================
#  Output
# =====================================================================
function write_csv(runs, points)
    mkpath(OUTDIR)
    open(joinpath(OUTDIR, "atlas_runs.csv"), "w") do io
        println(io, "point,kappa1,kappa2,C2,C3,ic,icarg,kind,mode,share,peaks,amp,corr_Crho,corr_CR,resid,t")
        for r in runs
            @printf(io, "%d,%d,%d,%.8g,%.8g,%s,%d,%s,%d,%.4g,%d,%.4g,%.4g,%.4g,%.3g,%.1f\n",
                    r.point, r.window..., r.C2, r.C3, r.ic, r.icarg, r.kind, r.nd, r.share,
                    r.peaks, r.amp, r.corr_Crho, r.corr_CR, r.resid, r.t)
        end
    end
    open(joinpath(OUTDIR, "atlas_points.csv"), "w") do io
        println(io, "point,kappa1,kappa2,C2,C3,linear_ok,nstar,growth_star,noise_modes,all_modes,nclasses,nshapes,any_unconverged")
        for p in points
            @printf(io, "%d,%d,%d,%.8g,%.8g,%s,%d,%.4g,%s,%s,%d,%d,%s\n",
                    p.point, p.window..., p.C2, p.C3, p.linear_ok, p.nstar, p.growth_star,
                    join(p.noise_modes, ";"), join(p.all_modes, ";"), p.nclasses, p.nshapes, p.any_unconverged)
        end
    end
end

const MODE_COLS = Dict(0 => :gray60, 1 => "#1b9e77", 2 => "#d95f02", 3 => "#7570b3", 4 => "#e7298a", 5 => "#66a61e")
modecol(n) = get(MODE_COLS, n, :black)

function region_outlines!(plt, R)
    for (w, P) in sort(collect(R))
        Q = vcat(P, P[1:1])
        plot!(plt, first.(Q), last.(Q); color = :black, lw = 0.8, label = "")
        c = (sum(first, P) / length(P), sum(last, P) / length(P))
        annotate!(plt, c[1], c[2], text("$(w[1])..$(w[2])", 7, :black))
    end
    plt
end

"Mode most often reached from noise at a point (ties -> smallest)."
function majority(v)
    isempty(v) && return 0
    u = sort(unique(v));  u[argmax([count(==(k), v) for k in u])]
end

function plot_atlas(points, R)
    mkpath(OUTDIR)
    base() = region_outlines!(plot(; xlabel = "C₂", ylabel = "C₃", legend = :outerright,
                                    size = (900, 650), titlefontsize = 10), R)
    # linear prediction
    pl = base(); title!(pl, "Linear prediction: fastest-growing mode n*")
    for n in sort(unique(p.nstar for p in points))
        ps = filter(p -> p.nstar == n, points)
        scatter!(pl, [p.C2 for p in ps], [p.C3 for p in ps]; color = modecol(n), ms = 4, msw = 0, label = "n* = $n")
    end
    # realised pattern from noise
    pr = base(); title!(pr, "Realised pattern from noise (majority class); ◆ = >1 class reached, × = unconverged")
    for n in sort(unique(majority(p.noise_modes) for p in points))
        ps = filter(p -> majority(p.noise_modes) == n, points)
        mono  = filter(p -> p.nclasses == 1, ps);  multi = filter(p -> p.nclasses > 1, ps)
        scatter!(pr, [p.C2 for p in mono], [p.C3 for p in mono]; color = modecol(n), ms = 4, msw = 0,
                 label = n == 0 ? "other uniform state" : "mode $n")
        scatter!(pr, [p.C2 for p in multi], [p.C3 for p in multi]; color = modecol(n), ms = 6, shape = :diamond,
                 msw = 0.5, label = "")
    end
    un = filter(p -> p.any_unconverged, points)
    isempty(un) || scatter!(pr, [p.C2 for p in un], [p.C3 for p in un]; shape = :x, color = :black, ms = 5, label = "unconverged")
    savefig(pl, joinpath(OUTDIR, "atlas_linear.png"))
    savefig(pr, joinpath(OUTDIR, "atlas_realised.png"))
    return pl, pr
end

"Per window: fraction of runs (noise ICs, cosine-seeded ICs) ending in each mode."
function plot_window_bars(runs)
    ws = sort(unique(r.window for r in runs))
    modes = sort(unique(r.nd for r in runs))
    panels = map((:perturb, :mode)) do ic
        M = [count(r -> r.window == w && r.ic === ic && r.nd == n, runs) /
             max(1, count(r -> r.window == w && r.ic === ic, runs)) for n in modes, w in ws]
        h = heatmap(1:length(ws), 1:length(modes), M; clims = (0, 1), color = :Blues,
                    xticks = (1:length(ws), ["$(w[1])..$(w[2])" for w in ws]),
                    yticks = (1:length(modes), [n == 0 ? "other uniform" : "mode $n" for n in modes]),
                    xlabel = "theorem window", title = ic === :perturb ? "noise ICs" : "cosine-seeded ICs")
        for i in eachindex(ws), k in eachindex(modes)
            M[k, i] > 0 && annotate!(h, i, k, text(@sprintf("%.2f", M[k, i]), 7, M[k, i] > 0.5 ? :white : :black))
        end
        h
    end
    plt = plot(panels...; layout = (2, 1), size = (900, 700))
    savefig(plt, joinpath(OUTDIR, "atlas_window_bars.png"))
    plt
end

"Final C, R, rho at each window's centroid, one panel per distinct state reached (up to reflection)."
function plot_examples(runs, finals, x)
    N = length(x)
    panels = []
    for w in sort(unique(r.window for r in runs))
        i = minimum(r.point for r in runs if r.window == w)          # centroid is the first point
        js = findall(r -> r.point == i, runs)
        reps = Int[]
        for j in js
            v = finals[j];  sc = max(maximum(abs, v), 1e-12)
            any(min(maximum(abs, v .- finals[k]), maximum(abs, reverse_fields(v, N) .- finals[k])) < 1e-2sc for k in reps) ||
                push!(reps, j)
        end
        for j in reps
            v = finals[j];  r = runs[j]
            ttl = @sprintf("%d..%d  %s → %s", w..., r.ic === :perturb ? "noise $(r.icarg)" : "cos $(r.icarg)x",
                           r.nd == 0 ? "other uniform" : "mode $(r.nd)")
            plt = plot(x, v[1:N]; color = CCOL, lw = 2, label = "Cdc42", title = ttl, titlefontsize = 8,
                       xlims = (0, π), xticks = ([0, π/2, π], ["0", "π/2", "π"]), legend = isempty(panels) ? :topright : false)
            plot!(plt, x, v[N+1:2N]; color = RCOL, lw = 2, label = "Rac")
            plot!(plt, x, v[2N+1:3N]; color = RHOCOL, lw = 2, label = "Rho")
            push!(panels, plt)
        end
    end
    plt = plot(panels...; layout = length(panels), size = (1400, 220 * cld(length(panels), 4)))
    savefig(plt, joinpath(OUTDIR, "atlas_examples.png"))
    plt
end
reverse_fields(v, N) = vcat(reverse(v[1:N]), reverse(v[N+1:2N]), reverse(v[2N+1:3N]))

# =====================================================================
if abspath(PROGRAM_FILE) == @__FILE__
    runs, points, R, finals, x = run_atlas()
    write_csv(runs, points)
    plot_atlas(points, R)
    plot_window_bars(runs)
    plot_examples(runs, finals, x)
    for w in sort(unique(p.window for p in points))
        ps = filter(p -> p.window == w, points)
        @printf("window %d..%d: n* = realised majority at %2d / %d points; >1 class at %2d\n", w...,
                count(p -> p.nstar == majority(p.noise_modes), ps), length(ps), count(p -> p.nclasses > 1, ps))
    end
    nbad = count(p -> !p.linear_ok, points)
    @printf("linear window mismatches: %d / %d points\n", nbad, length(points))
    println("figures and CSVs in ", OUTDIR)
end
