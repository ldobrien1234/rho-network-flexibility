# =====================================================================
#  run_simulations.jl
#
#  Demonstrates the central claim of the companion LaTeX section:
#
#   * REGIME (P)  Rac -> Rho  (activating,  f_rhoR > 0):  the reduced
#     Cdc42-Rac-Rho network is a FLEXIBLE pattern former.  The same
#     interaction graph, with the same fixed signs, can be tuned so that
#     ANY prescribed admissible wavenumber n = 1,2,3,4 is the primary
#     (Turing) instability.  Four parameter sets, built by the backward
#     construction of Proposition "Existence of a primary instability at
#     a prescribed mode", each select a different number of stripes.
#
#   * REGIME (S)  Rac -| Rho  (inhibiting,  f_rhoR < 0):  the network is
#     INFLEXIBLE.  By Theorem "Patterning forces C_{Crho} > 0 and C_3 < 0"
#     and its Corollary, NO parameter choice yields a finite-wavelength
#     instability.  Weak feedback -> a stable homogeneous state
#     (perturbations decay); strong feedback -> a bistable, winner-take-all
#     SWITCH (uniform final state, position/identity set by initial data).
#     Never a periodic pattern.
#
#   * A random parameter scan (scan_regime) makes the dichotomy
#     quantitative: a Turing instability is found for a healthy fraction
#     of stable-homogeneous regime-(P) draws and for exactly zero
#     regime-(S) draws.
#
#  Output: PNG figures + a CSV summary in ./output/.
# =====================================================================

ENV["GKSwstype"] = "100"          # headless GR
using Plots
using Printf
using DelimitedFiles
include("cdc42_rac_rho.jl")

import Plots: mm
gr(size = (900, 600), dpi = 130, legend = :outerright, framestyle = :box,
   titlefontsize = 10, guidefontsize = 9, legendfontsize = 8,
   left_margin = 6mm, bottom_margin = 6mm, top_margin = 3mm)
const OUT = joinpath(@__DIR__, "output")
mkpath(OUT)

# palette
const CCOL = "#1b7837"; const RCOL = "#2166ac"; const RHOCOL = "#b2182b"

# ---------------------------------------------------------------------
#  helpers
# ---------------------------------------------------------------------
interior_extrema(v) = count(i -> (v[i]-v[i-1])*(v[i+1]-v[i]) < 0, 2:length(v)-1)

# dominant admissible wavenumber of a profile on [0,2pi] (no-flux) via cosine content
function dominant_mode(v; nmax = 20)
    w = v .- sum(v)/length(v)
    x = range(0, 2π; length = length(v))
    best_n, best_p = 0, 0.0
    for n in 1:nmax
        c = sum(w .* cos.(n .* x))
        c^2 > best_p && (best_p = c^2; best_n = n)
    end
    best_n
end

function disp_panel(d; title = "", target = -1, nshow = 9, zoom = true)
    sel  = d.n .<= nshow
    ns, gs = d.n[sel], d.growth[sel]
    bar_c = [g > 1e-7 ? :crimson : :steelblue for g in gs]
    gmax = maximum(d.growth)
    yl = zoom ? (-0.05, max(0.008, 3gmax)) : (min(1.15*minimum(gs), -0.01), max(0.02, 1.15gmax))
    plt = bar(ns, gs; label = "", color = bar_c, bar_width = 0.62,
              xlabel = "wavenumber n  (mode cos nx)", ylabel = "max Re λ(n)",
              title = title, legend = false, ylims = yl)
    hline!(plt, [0.0]; color = :black, lw = 1.2, label = "")
    target > 0 && vline!(plt, [target]; color = :seagreen, ls = :dash, lw = 2, label = "")
    plt
end

function profile_panel(r; title = "", showlegend = true, ylims = :auto)
    plt = plot(r.x, r.C;   color = CCOL,   lw = 2.4, label = "Cdc42",
               xlabel = "x", ylabel = "activity", title = title, ylims = ylims,
               legend = showlegend ? :outerright : false, xlims = (0, 2π),
               xticks = ([0, π, 2π], ["0", "π", "2π"]))
    plot!(plt, r.x, r.R;   color = RCOL,   lw = 2.4, label = "Rac")
    plot!(plt, r.x, r.rho; color = RHOCOL, lw = 2.4, label = "Rho")
    # homogeneous state reference
    hline!(plt, [r.ustar[1]]; color = CCOL,   ls = :dot, lw = 1, label = "")
    hline!(plt, [r.ustar[2]]; color = RCOL,   ls = :dot, lw = 1, label = "")
    hline!(plt, [r.ustar[3]]; color = RHOCOL, ls = :dot, lw = 1, label = "")
    plt
end

# accumulate a summary table
rows = Vector{Vector{Any}}()
push_row!(tag, p, d, r, note) = push!(rows, Any[
    tag, p.regime, p.target_mode,
    round(p.dC; sigdigits = 3), round(p.dR; sigdigits = 3),
    round(d.growth0; digits = 4), d.homog_stable, d.nstar,
    round(d.growth_star; digits = 4), d.turing,
    r === nothing ? "" : interior_extrema(r.rho),
    r === nothing ? "" : round(minimum(r.rho); digits = 3),
    r === nothing ? "" : round(maximum(r.rho); digits = 3),
    r === nothing ? "" : round(r.resid; sigdigits = 2), note])

# =====================================================================
#  1.  REGIME (P):  flexibility -- select modes n = 1,2,3,4
# =====================================================================
println("="^70, "\n REGIME (P):  Rac -> Rho activating  --  flexible pattern former\n", "="^70)

P_targets = 1:5          # constructed + linear dispersion for all of these
# full PDE integration: target mode => (grid N, integration time).  High n needs
# a long transient because near onset the growth rate is O(1e-3) (the high-n well
# is shallow: |g''(-n^2)| <~ 18/n^4 while the homogeneous state stays stable).
P_sim_plan = Dict(1 => (161, 6000.0), 2 => (161, 12000.0), 3 => (161, 55000.0))
const P_GAIN = 1.004     # small, fixed supercriticality (weakly nonlinear -> mode stays put)
P_params  = RDParams[]
P_disp    = []
P_sims    = Dict{Int,Any}()

for n in P_targets
    p0, info = construct_regime_P(n)
    pon = apply_twocycle_gain(p0, P_GAIN)
    d   = dispersion(pon; nmax = 16)
    nun = count(>(1e-7), d.growth[2:end])
    @printf("n*=%d :  synthesised dC=drho=%.4g, dR=%.1f ; 2-cycle C_Crho=%.3f (>0), 3-cycle C_3=%.3f (<0)\n",
            n, info.dC, pon.dR, info.CCrho, info.C3)
    @printf("        supercritical (2-cycle gain %.3f): homog_stable=%s, fastest mode n=%d, growth=%.5f  [# unstable modes = %d]\n",
            P_GAIN, d.homog_stable, d.nstar, d.growth_star, nun)
    push!(P_params, pon); push!(P_disp, d)
    if haskey(P_sim_plan, n)
        N, tmx = P_sim_plan[n]
        r = simulate(pon; N = N, tmax = tmx, dt = 5e-3, seed = 20 + n,
                     amp = 1e-3, tol = 1e-8)
        P_sims[n] = r
        dm = dominant_mode(r.rho)
        ok = dm == n ? "OK" : "!! got $dm"
        @printf("        steady state (random IC): dominant wavenumber = %d (target %d) [%s], Rho in [%.3f, %.3f], resid=%.1e, t=%.0f\n\n",
                dm, n, ok, minimum(r.rho), maximum(r.rho), r.resid, r.t)
        push_row!("P$n", pon, d, r, "primary Turing mode n=$n (random IC -> stationary mode $dm)")
        fig = plot(disp_panel(d; title = "dispersion  (regime P, tuned for n=$n)", target = n),
                   profile_panel(r; title = "steady state  (stationary wavenumber $dm)");
                   layout = (1, 2), size = (1150, 430))
        savefig(fig, joinpath(OUT, "P$(n)_mode$(n).png"))
    else
        println()
        tag = nun == 1 ? "unique unstable admissible mode n=$n" : "fastest-growing mode n=$n"
        push_row!("P$n", pon, d, nothing, "linear analysis only: $tag")
    end
end
P_sim_targets = sort(collect(keys(P_sims)))

# combined dispersion figure: the peak walks across n  (full + zoom on the crossing)
function disp_family(ylims, title, leg)
    p = plot(; xlabel = "admissible wavenumber n", ylabel = "max Re λ(n)  (growth rate)",
             title = title, legend = leg, xlims = (-0.3, 8.4), ylims = ylims, xticks = 0:8)
    for (k, n) in enumerate(P_targets)
        d = P_disp[k]; sel = d.n .<= 8
        plot!(p, d.n[sel], d.growth[sel]; lw = 2, marker = :circle, ms = 4,
              label = "tuned for n = $n")
        gk = d.growth[n + 1]
        (ylims[1] <= gk <= ylims[2]) && scatter!(p, [n], [gk]; ms = 9, msw = 2,
              markershape = :star5, color = :gold, label = "")
    end
    hline!(p, [0.0]; color = :black, lw = 1.5, label = "")
    p
end
figD = plot(disp_family((-0.45, 0.02), "full range", false),
            disp_family((-0.03, 0.006), "zoom on the zero crossing  (★ = target mode)", :bottomleft);
            layout = (1, 2), size = (1300, 500),
            plot_title = "Regime (P): one graph, one set of signs — the dispersion peak is tunable to any n",
            plot_titlefontsize = 12)
savefig(figD, joinpath(OUT, "P_dispersion_family.png"))

# montage of the steady states
nsim = length(P_sim_targets)
panelsM = map(enumerate(P_sim_targets)) do (k, n)
    r = P_sims[n]
    p = plot(r.x, r.C; color = CCOL, lw = 2, label = k == 1 ? "Cdc42" : "",
             title = "tuned for n = $n  (stationary wavenumber $(dominant_mode(r.rho)))",
             titlefontsize = 9, xlims = (0, 2π), xticks = ([0, π, 2π], ["0", "π", "2π"]),
             xlabel = "x", ylabel = k == 1 ? "activity" : "",
             legend = k == 1 ? :top : false)
    plot!(p, r.x, r.R;   color = RCOL,   lw = 2, label = k == 1 ? "Rac" : "")
    plot!(p, r.x, r.rho; color = RHOCOL, lw = 2, label = k == 1 ? "Rho" : "")
    p
end
figM = plot(panelsM...; layout = (1, nsim), size = (450nsim, 420), top_margin = 9mm,
            plot_title = "Regime (P): one graph, one set of signs — flexible wavelength selection on [0, 2π]",
            plot_titlefontsize = 12)
savefig(figM, joinpath(OUT, "P_steady_states_montage.png"))

# =====================================================================
#  2.  REGIME (S):  inflexibility
# =====================================================================
println("="^70, "\n REGIME (S):  Rac -| Rho inhibiting  --  inflexible (switch, no pattern)\n", "="^70)

# 2a. mild feedback -> stable homogeneous state, all modes decay
pS_stable = build_regime_S(; strength = 0.6, label = "S (weak): stable homogeneous")
dS_stable = dispersion(pS_stable; nmax = 16)
rS_stable = simulate(pS_stable; N = 161, tmax = 4000.0, dt = 5e-3, seed = 1, amp = 0.08, tol = 1e-7)
@printf("weak (strength 0.6):  homog_stable=%s, max finite-wavelength growth=%.4f (<0), turing=%s\n",
        dS_stable.homog_stable, maximum(dS_stable.growth[2:end]), dS_stable.turing)
@printf("   -> simulation returns to homogeneous:  Rho in [%.4f, %.4f]  (u* = %.3f)\n\n",
        minimum(rS_stable.rho), maximum(rS_stable.rho), rS_stable.ustar[3])
push_row!("S_stable", pS_stable, dS_stable, rS_stable, "weak feedback: perturbations decay")

# 2b. strong feedback -> bistable winner-take-all switch (all loops positive => multistable)
pS_switch = build_regime_S(; strength = 2.0, m = 6, KC = 0.8, Krho1 = 0.8, Krho2 = 0.8,
                             KR = 6.0, label = "S (strong): winner-take-all switch")
dS_switch = dispersion(pS_switch; nmax = 16)
cs = count_homog_states(pS_switch)
rS_A    = simulate(pS_switch; N = 161, tmax = 4000.0, dt = 3e-3, seed = 1, ic = :perturb, amp = 0.03, tol = 1e-7)
rS_B    = simulate(pS_switch; N = 161, tmax = 4000.0, dt = 3e-3, seed = 4, ic = :perturb, amp = 0.03, tol = 1e-7)
rS_step = simulate(pS_switch; N = 161, tmax = 4000.0, dt = 3e-3, seed = 1, ic = :step,    tol = 1e-7)
@printf("strong (strength 2.0, m 6):  turing=%s ; %d homogeneous states, %d stable:\n",
        dS_switch.turing, cs.nstates, cs.nstable)
for s in cs.states; @printf("       (C,R,Rho) = %s\n", string(round.(s; digits = 3))); end
@printf("   IC noise seed 1 -> uniform Rho = %.3f ;  seed 4 -> uniform Rho = %.3f ;  step IC -> uniform Rho = %.3f\n",
        sum(rS_A.rho)/length(rS_A.rho), sum(rS_B.rho)/length(rS_B.rho), sum(rS_step.rho)/length(rS_step.rho))
println("   => winner-take-all: final state is UNIFORM and initial-condition dependent; no wavelength.\n")
push_row!("S_switch_seedA", pS_switch, dS_switch, rS_A, "bistable switch: noise IC -> uniform winner A")
push_row!("S_switch_seedB", pS_switch, dS_switch, rS_B, "bistable switch: noise IC -> uniform winner B")

sA = disp_panel(dS_stable; title = "weak feedback: every admissible mode decays (max Re λ < 0)", zoom = false)
sB = profile_panel(rS_stable; title = "weak S: perturbation decays back to homogeneous",
                   showlegend = true, ylims = (0, 1.6))
sC = disp_panel(dS_switch; title = "strong feedback: unstable, but the peak is at n → 0, not finite n")
sD = plot(rS_A.x, rS_A.rho; color = RHOCOL, lw = 2.5, label = "Rho — noise seed 1",
          xlabel = "x", ylabel = "Rho activity", xlims = (0, 2π), ylims = (-0.2, 3.9),
          xticks = ([0, π, 2π], ["0", "π", "2π"]), legend = :right,
          title = "strong S: winner-take-all — uniform, IC-dependent")
plot!(sD, rS_B.x, rS_B.rho; color = :darkorange, lw = 2.5, ls = :dash, label = "Rho — noise seed 4")
figS = plot(sA, sB, sC, sD; layout = (2, 2), size = (1250, 820),
            plot_title = "Regime (S): switch or nothing — never a finite-wavelength pattern",
            plot_titlefontsize = 12)
savefig(figS, joinpath(OUT, "S_switch_and_stable.png"))

# side-by-side: flip the single Rac->Rho edge of a working P set  (plot Rho only)
p_flip = flip_to_regime_S(P_params[1]; note = "P1 with Rac->Rho edge negated")
d_flip = dispersion(p_flip; nmax = 16)
r_flip = simulate(p_flip; N = 161, tmax = 3000.0, dt = 5e-3, seed = 5, amp = 0.05, tol = 1e-7)
push_row!("S_flip_of_P1", p_flip, d_flip, r_flip, "negate one edge: pattern -> switch")
r1 = P_sims[1]
yl = (0, 1.15 * max(maximum(r1.rho), maximum(r_flip.rho)))
fL = plot(r1.x, r1.rho; color = RHOCOL, lw = 2.8, label = "", xlims = (0, 2π), ylims = yl,
          xticks = ([0, π, 2π], ["0", "π", "2π"]), xlabel = "x", ylabel = "Rho activity",
          title = "P1: Rac → Rho activating\nstationary Turing pattern (n = 1)")
fR = plot(r_flip.x, r_flip.rho; color = RHOCOL, lw = 2.8, label = "", xlims = (0, 2π), ylims = yl,
          xticks = ([0, π, 2π], ["0", "π", "2π"]), xlabel = "x", ylabel = "Rho activity",
          title = "same rate constants, Rac ⊣ Rho negated\nuniform (winner-take-all), no pattern")
figF = plot(fL, fR; layout = (1, 2), size = (1200, 440),
            plot_title = "Negate the single Rac→Rho edge: the finite-wavelength pattern is gone",
            plot_titlefontsize = 12)
savefig(figF, joinpath(OUT, "edge_flip_comparison.png"))

# =====================================================================
#  3.  Structural scan:  the dichotomy is not about tuning
# =====================================================================
println("="^70, "\n STRUCTURAL SCAN (random sign-legal parameters)\n", "="^70)
sP = scan_regime(:P; ndraw = 6000, seed = 11)
sS = scan_regime(:S; ndraw = 6000, seed = 11)
for s in (sP, sS)
    mx = isempty(s.maxg_turing) ? NaN : maximum(s.maxg_turing)
    @printf("regime %s :  %d draws, %d with stable homogeneous state, %d of those Turing-unstable (%.2f%%);  max finite-wavelength growth = %+.3f\n",
            s.regime, s.ndraw, s.nstable, s.nturing, 100s.frac_turing, mx)
end

bins = range(-1.05, 0.65; length = 44)
figH = histogram(sP.maxg_turing; bins = bins, alpha = 0.55,
                 label = "regime P  (Rac → Rho, +) :  $(sP.nturing) / $(sP.nstable) are Turing-unstable",
                 color = :seagreen,
                 xlabel = "max finite-wavelength growth rate   maxₙ≥₁ Re λ(n)   (over stable-homogeneous draws)",
                 ylabel = "count", legend = :topright,
                 title = "Structural dichotomy: with Rac ⊣ Rho, no random draw ever crosses zero")
histogram!(figH, sS.maxg_turing; bins = bins, alpha = 0.55, color = :firebrick,
           label = "regime S  (Rac ⊣ Rho, −) :  $(sS.nturing) / $(sS.nstable) are Turing-unstable")
vline!(figH, [0.0]; color = :black, lw = 2.5, ls = :dash, label = "Turing threshold")
# headroom so the tallest bar is not clipped
_hc = maximum(step -> count(g -> bins[step] <= g < bins[step+1], sP.maxg_turing), 1:length(bins)-1)
plot!(figH; size = (1050, 640), ylims = (0, 1.1_hc))
savefig(figH, joinpath(OUT, "structural_scan_histogram.png"))

# =====================================================================
#  4.  write summary CSV + console table
# =====================================================================
header = ["set" "regime" "target_n" "dC" "dR" "growth0" "homog_stable" "n_star" "growth_star" "turing" "rho_interior_extrema" "rho_min" "rho_max" "residual" "note"]
open(joinpath(OUT, "summary.csv"), "w") do io
    writedlm(io, vcat(header, permutedims(hcat(rows...))), ',')
end

println("\n", "="^70, "\n SUMMARY\n", "="^70)
println(rpad("set", 16), rpad("regime", 8), rpad("target_n", 9), rpad("homog_stbl", 12),
        rpad("n*", 5), rpad("growth*", 11), rpad("turing", 8), rpad("Rho extrema", 12), "note")
for row in rows
    println(rpad(row[1], 16), rpad(string(row[2]), 8), rpad(string(row[3]), 9),
            rpad(string(row[7]), 12), rpad(string(row[8]), 5),
            rpad(string(row[9]), 11), rpad(string(row[10]), 8),
            rpad(string(row[11]), 12), row[15])
end

println("\nFigures written to: ", OUT)
foreach(f -> println("  ", f), sort(readdir(OUT)))
