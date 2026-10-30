# Flexibility vs. inflexibility of patterning in the Cdc42–Rac–Rho network

A standard three-node reaction–diffusion model is built. We simulate theoretical predictions of pattern formation, which are based on a model-independent and network based approach.


Julia 1.10+ (tested on 1.12).  Dependencies: `Plots`, `ForwardDiff` (both in the
default environment here).  Run:

```
julia run_simulations.jl            # figures + summary.csv in output/
julia -t auto cycle_atlas.jl        # attractor atlas, output/atlas/
julia -t 3 onset_sweeps.jl          # onset and Hopf sweeps, output/sweeps/
```

`plot_parameter_region.m` (MATLAB) draws the theorem regions in `(C2, C3)`;
edit the parameter block at the top and run it.

---

## The model  (`cdc42_rac_rho.jl`)

`u = (C, R, ρ)` = active Cdc42, Rac, Rho on `x ∈ [0, π]`, no-flux:

```
∂_t u = D ∂_xx u + F(u)

f_C = b_C   + a_C   · hillI(ρ, K_C ,m)                              − γ_C  C
f_R = b_R   + a_R   · hillA(C, K_R ,m)                              − γ_R  R
f_ρ = b_ρ   + a_ρ1  · hillI(C, K_ρ1,m) + s · a_ρ2 · hillX(R,K_ρ2,m) − γ_ρ  ρ
```

`hillA` / `hillI` are increasing / decreasing Hill functions.  Each node has
linear turnover (`f_ii = −γ_i < 0`) and additive Hill regulation — the minimal
smooth realisation consistent with the interaction graph.  The interaction signs
are exactly those of the text:

| edge | sign | mechanism |
|---|---|---|
| Cdc42 → Rac | `f_RC > 0` | Cdc42 activates Rac |
| Rho ⊣ Cdc42, Cdc42 ⊣ Rho | `f_Cρ, f_ρC < 0` | mutual inhibition |
| **Rac → Rho** | **`f_ρR` : sign is the control parameter** | see below |

**Regime P** (`hillX = hillA`, `f_ρR > 0`): Rac *activates* Rho — Rac recruits
Arhgef11/Arhgef12 (Lbc-type GEFs) as effectors (Marston et al., *Nat Commun*
2023).

**Regime S** (`hillX = hillI`, `f_ρR < 0`): Rac *inhibits* Rho — Rac → Pak1/Pak4
phosphorylates and inactivates RhoA GEFs p115-RhoGEF / GEF-H1 / PDZ-RhoGEF
(Rosenfeldt et al. 2006; Guilluy et al. 2011).

Dispersion: admissible modes are `cos(n x)`, `n = 0,1,2,…` (`B = −n²`), and the
growth rate of mode `n` is `maxRe eig(J − n² D)` with `J = DF(u*)`.

## What the scripts show

### `construct_regime_P(n)`  — the paper's backward construction

Given a target mode `n`, fixes the *root geometry* of the dispersion cubic
`g(B) = det L(B)` so that `B = −n²` is a **double root** (marginal), solves
backward for the cycle quantities `C_{Cρ} = f_Cρ f_ρC` (2-cycle) and
`C_3 = f_RC f_ρR f_Cρ` (3-cycle), and realises them with sign-legal Hill
kinetics at `u* = (1,1,1)`.  Output confirms `C_{Cρ} > 0` and `C_3 < 0` for every
`n` (Theorem "Patterning forces `C_{Cρ} > 0` and `C_3 < 0`").

`push_to_onset` then raises a single 2-cycle gain infinitesimally past the
marginal point to get a genuine one-mode instability.

### Regime P is flexible

`P_dispersion_family.png` — five constructed sets; the linear dispersion peak
sits **exactly** on `n = 1,2,3,4,5` in turn (right panel, ★).  Same graph, same
signs, only rate constants change.  This is the paper's actual claim — a
statement about the *primary instability*, i.e. a linear one.

`P_steady_states_montage.png`, `P{1,2,3}_mode{1,2,3}.png` — PDE integration from
small random noise, at a fixed small supercriticality (2-cycle gain 1.004),
converges to a **stationary `n`-stripe pattern** (dominant cosine wavenumber
`= n`, verified by `dominant_mode`).

`n = 4, 5` are shown by the dispersion relation only.  With a linearly stable
homogeneous state the well curvature at mode `n` obeys `|g''(−n²)| ≲ 18/n⁴`, so
the high-`n` wells are very flat: nonlinear selection there needs the capture
band to be tuned narrower than the lattice gap (the text's onset proposition),
and away from that narrow window a Turing system on a finite interval coarsens
toward the longest-wavelength large-amplitude state.  `n = 1,2,3` are safely
inside the robust window; `n = 4,5` are not, in this three-node model with an
order-one homogeneous decay rate.

### Regime S is inflexible

`S_switch_and_stable.png`:

* **weak feedback** → stable homogeneous state, **every** admissible mode decays
  (`turing = false`), a perturbation relaxes back to homogeneous;
* **strong feedback** → the homogeneous state loses stability, but the growth
  rate is *maximal at `n → 0`* (spatially uniform), not at any finite `n`.  The
  network is bistable and goes **winner-take-all**: the final state is uniform
  and depends on the initial condition (`S_switch_and_stable.png`, bottom-right:
  two noise seeds → two different uniform Rho levels).

`edge_flip_comparison.png` — take a working regime-P set and **negate the single
Rac → Rho edge**: the `n = 1` Turing pattern is replaced by a uniform
winner-take-all state.

### The dichotomy is structural, not fine-tuning

`structural_scan_histogram.png` / `scan_regime` — 6000 random sign-legal
parameter draws per regime.  Among draws with a linearly stable homogeneous
state:

```
regime P :  ~6 %   admit a finite-wavelength (Turing) instability;  max growth ≈ +0.5
regime S :   0     (0 / ~1600);                                     max growth < 0  (bounded away from 0)
```

This is the numerical form of the Corollary (Regime dichotomy): with Rac ⊣ Rho
the 3-cycle `C_3 > 0` is forced, the necessary condition `C_3 < 0` can never be
met, and no amount of parameter tuning yields a pattern.

## Theorem regions in cycle space

The theorems are stated in the cycle quantities `C2 = f_Cρ f_ρC` and
`C3 = f_RC f_ρR f_Cρ`, with the diagonal `f_ii` and diffusivities `d_i > 0` held
fixed.  The operator is `J − k² D`, so mode `k` corresponds to `Λ = −k²`.

### `plot_parameter_region.m`

For a target window of unstable wavenumbers `k = κ1..κ2` (and every sub-window),
draws the Theorem A region in `(C2, C3)`.  The variables are

```
φ_i = f_ii / d_i,    S = −Σ φ_i,    ξ2 = Σ_{i<j} φ_i φ_j,    ξ3 = Π φ_i
```

(`f_ii` may be any real number; `f_ii > 0` lets `S < 0`, which reaches the empty,
triangle and first quadrilateral rows of the geometry table).  Each polygon vertex is

```
C2 = d_C d_ρ (ξ2 − x y − (x + y) z),     z = S − x − y
C3 = d_C d_R d_ρ (φ_R C2/(d_C d_ρ) − x y z − ξ3)
```

The `−x y` term is intended, and the region area jumps at `S = m + M` and
`S = σ + M` are real.  Theorem A regions are also intersected with the half-plane

```
C3 > −(f_CC + f_ρρ) C2 + (f_CC + f_RR)(f_CC + f_ρρ)(f_RR + f_ρρ)
```

and the window `κ1..κ2` is guaranteed only on **polygon ∩ half-plane**, not on
the polygon alone.  Its boundary `W₀` is the Routh–Hurwitz condition for `J`
itself; outside it the uniform mode `k = 0` goes unstable as a complex pair.
Theorem B (`θ_ij = (f_ii + f_jj)/(d_i + d_j)`) is kept in the file but commented
out; it has no half-plane condition.

Figure 1 shows one subplot per window, and Figure 2 overlays every region with
the line `L : C3 = d_R φ_R C2 − d_C d_R d_ρ ξ3`, which every region lies on or below.

### `cycle_atlas.jl` — what the nonlinear system does inside each region

Uses the same fixed diagonal as the MATLAB defaults,
`f = (−5, −1, −5)`, `d = (0.05, 0.5, 0.05)`, on `x ∈ [0, π]`, for every window
with `κ2 ≤ 4`.  It works in four steps:

1. Builds each region as polygon ∩ half-plane (clipped to `C2 ≥ 0`, `C3 ≤ 0`).
2. Samples 30 points per region, starting with the centroid.
3. Lifts each point to sign-legal Hill kinetics whose Jacobian at
   `u* = (1,1,1)` has exactly that `(C2, C3)` (`lift`).
4. Checks that the unstable modes of `J − k² D` are exactly `κ1..κ2`.  It then
   integrates the PDE from three random-noise ICs and from one `cos(k x)`-seeded
   IC per mode in the window, and classifies each final state as mode-`n`
   pattern, `u*`, another uniform state, or unconverged.

In the current output, the linear window matched at all 210 points.  The
single-mode windows `1..1` and `2..2` always reach one class of state.  Wide
windows are multistable: more than one class of state is reached at 26/30 points
of `1..2` and `1..3` and at 30/30 points of `1..4`.  So which pattern appears
depends on the initial condition, not only on the fastest linear mode.

| output | contents |
|---|---|
| `atlas_points.csv` | one row per point: window, `(C2, C3)`, `linear_ok`, fastest mode `n*`, modes reached from noise / from all ICs, number of distinct states |
| `atlas_runs.csv` | one row per simulation: IC, final class, dominant mode and its share, peak count, amplitude, correlations, residual |
| `atlas_linear.png` | sampled points coloured by linear fastest mode `n*` |
| `atlas_realised.png` | same points coloured by the majority mode reached from noise (◆ = more than one class reached) |
| `atlas_window_bars.png` | per window, fraction of runs ending in each mode (noise ICs vs cosine-seeded ICs) |
| `atlas_examples.png` | final `C, R, ρ` profiles at each window's centroid, one panel per distinct state |

### `onset_sweeps.jl` — crossing the boundaries slowly

**Turing onset.**  For `k = 1, 2, 3` the `k..k` region has an edge on the neutral
line of mode `k`,

```
T_k :  g_k(C2, C3) = Π_i (f_ii − d_i k²) − C2 (f_RR − d_R k²) + C3 = 0
```

The path starts at the midpoint of that edge and runs along its unit normal, so
`s` is the signed distance past onset: `s < 0` is stable, and for `s > 0` only
mode `k` is unstable.  The parameters are stepped forward and then backward.  At
each step the kinetics are re-lifted, and the previous profile plus a `10⁻⁶` kick
is relaxed to steady state.  A smooth (supercritical) onset grows like `√s` and
retraces itself.  A sudden (subcritical) onset jumps to finite amplitude and
persists below `s = 0` on the way back.

In the current output, all three crossings show hysteresis.  The pattern (or
the state it leads to) survives on the way back to `s ≈ −0.08` for `k = 1`,
`−0.23` for `k = 2` and `−0.49` for `k = 3`, against steps of about `0.02`.
Two other results:
* For `k = 1`, the mode-1 pattern appears just past onset, but further out it
  gives way to a uniform Rho-high state.
* For `k = 3`, the state reached is a **mode-2** pattern, even though only mode 3
  is linearly unstable.

**Hopf crossing.**  The second sweep crosses the half-plane boundary `W₀` at
`C2 ≈ 33`.  Outside the half-plane the network oscillates in time while staying
spatially uniform.  In the current output, the oscillation starts at the first
step past `W₀` and its amplitude grows like `√s` (fitted in the figure).  Its
period is close to the linear value `2π/ω`.  On the way back the oscillation
disappears within one step of `W₀`, so the Hopf onset looks supercritical.

| output | contents |
|---|---|
| `sweeps.csv` | every step of the three Turing sweeps: `s`, `(C2, C3)`, growth of mode `k`, state, dominant mode, means, residual |
| `onset_sweeps.png` | pattern amplitude vs `s`, forward ● / backward ○, coloured by state |
| `onset_profiles.png` | profiles at first step past onset, end of the forward sweep, and last pattern on the way back |
| `hopf_sweep.csv`, `hopf_sweep.png` | Hopf sweep: amplitude and period vs `s`, with time series of the spatial means |

## Files

| file | contents |
|---|---|
| `cdc42_rac_rho.jl` | model, Jacobian, homogeneous solver, linear dispersion, IMEX PDE solver, `construct_regime_P`, `push_to_onset`, `flip_to_regime_S`, `build_regime_S`, `scan_regime` |
| `run_simulations.jl` | builds all parameter sets, simulates, writes figures + `summary.csv` |
| `cycle_atlas.jl` | Theorem A regions (polygon ∩ half-plane), `lift` from `(C2, C3)` to Hill kinetics, attractor atlas |
| `onset_sweeps.jl` | parameter sweeps across `T_k` (Turing) and `W₀` (Hopf); includes `cycle_atlas.jl` |
| `plot_parameter_region.m` | MATLAB plots of the Theorem A regions and half-plane (Theorem B commented out) |
| `output/summary.csv` | one row per parameter set: regime, homogeneous stability, fastest mode `n*`, its growth rate, `turing` flag, realised pattern wavenumber |
