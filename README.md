# Flexibility vs. inflexibility of patterning in the Cdc42–Rac–Rho network

Numerical companion to the LaTeX section *"Flexibility of Spatial Patterning in
the Cdc42–Rac–Rho Network."*  A standard three-node reaction–diffusion model is
built to match the text, and parameter sets are found that show

* the **flexible** ("patterning", regime **P**) network selecting several
  different wavelengths, and
* the **inflexible** ("switch", regime **S**) network never producing a
  finite-wavelength pattern, for any parameters.

Julia 1.10+ (tested on 1.12).  Dependencies: `Plots`, `ForwardDiff` (both in the
default environment here).  Run:

```
julia run_simulations.jl
```

Figures and `summary.csv` are written to `output/`.

---

## The model  (`cdc42_rac_rho.jl`)

`u = (C, R, ρ)` = active Cdc42, Rac, Rho on `x ∈ [0, 2π]`, no-flux:

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

## Files

| file | contents |
|---|---|
| `cdc42_rac_rho.jl` | model, Jacobian, homogeneous solver, linear dispersion, IMEX PDE solver, `construct_regime_P`, `push_to_onset`, `flip_to_regime_S`, `build_regime_S`, `scan_regime` |
| `run_simulations.jl` | builds all parameter sets, simulates, writes figures + `summary.csv` |
| `output/summary.csv` | one row per parameter set: regime, homogeneous stability, fastest mode `n*`, its growth rate, `turing` flag, realised pattern wavenumber |
