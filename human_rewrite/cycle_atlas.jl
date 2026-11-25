# =====================================================================
#  cycle_atlas.jl
#
#  Checks the theoretical partition of cycle space (C_2, C_3), where
#
#    C_2 = f_Crho f_rhoC           (Cdc42 <-> Rho 2-cycle)
#    C_3 = f_RC f_rhoR f_Crho      (Cdc42 -> Rac -> Rho -| Cdc42 3-cycle)
#
#  For each band of modes, the theory predicts a region of (C_2, C_3)
#  in which exactly those modes are unstable. For points sampled inside
#  each region, we lift to parameters for the model in cdc42_rac_rho.jl.
#  Then we check
#
#    (1) the unstable modes at the homogeneous steady-state.
#    (2) the steady states from different initial conditions
# =====================================================================

include("cdc42_rac_rho.jl")

using LazySets
using Random
using Plots
using ForwardDiff
using LinearAlgebra

## -------------------------- fixed parameters ------------------------
# Diagonal of the Jacobian
const f_CC     = -g_C
const f_RR     = -g_R
const f_rhorho = -g_rho

# phi_i = f_ii / d_i and S_phi = -(phi_C + phi_R + phi_rho)
const phi_C   = f_CC     / d_C
const phi_R   = f_RR     / d_R
const phi_rho = f_rhorho / d_rho
const S_phi   = -(phi_C + phi_R + phi_rho)

# theta_ij = (f_ii + f_jj) / (d_i + d_j) and S_theta = -(theta_CR + theta_Crho + theta_Rrho)
const theta_CR   = (f_CC + f_RR)     / (d_C + d_R)
const theta_Crho = (f_CC + f_rhorho) / (d_C + d_rho)
const theta_Rrho = (f_RR + f_rhorho) / (d_R + d_rho)
const S_theta    = -(theta_CR + theta_Crho + theta_Rrho)

# Max and minimum wave numbers to consider
const κ_min = 1
const κ_max = 6

# Off-diagonal Jacobian entry Cdc42 -> Rac, held fixed in the lift
const f_RC = 1.0

## -------------------------- vertex pairs ----------------------------
"""
    band_modes(k_1, k_2)

Modes mu = -k^2 of the domain [0, pi] around the band k_1..k_2, in
increasing order m < sigma_1 <= sigma_n < M:

    m       = -(k_2 + 1)^2     first stable mode below the band
    sigma_1 = -k_2^2           lowest mode in the band
    sigma_n = -k_1^2           highest mode in the band
    M       = -(k_1 - 1)^2     first stable mode above the band (0 if k_1 = 1)
"""
function band_modes(k_1, k_2)
    m       = -(k_2 + 1)^2
    sigma_1 = -k_2^2
    sigma_n = -k_1^2
    M       = -(k_1 - 1)^2
    return m, sigma_1, sigma_n, M
end

"""
    vertex_pairs(S, m, sigma_1, sigma_n, M)

The pairs (x, y) of the vertices of the region for a band, read off from
the table by the value of S (S_phi for the steady-state theorem, S_theta
for the Hopf theorem).  Each row of the table is one branch below.  The thresholds are in increasing order because
m + M < sigma_1 + sigma_n for the modes -k^2.
"""
function vertex_pairs(S, m, sigma_1, sigma_n, M)
    if S <= m + sigma_n
        # empty
        return []
    elseif S < m + M
        # triangle
        return [(m, sigma_n), (m, 0), (sigma_n, 0)]
    elseif S < sigma_1 + sigma_n
        # quadrilateral
        return [(m, M), (m, sigma_n), (M, 0), (sigma_n, 0)]
    elseif S < sigma_1 + M
        # pentagon
        return [(m, M), (m, sigma_n), (sigma_1, sigma_n), (sigma_1, 0), (M, 0)]
    else
        # quadrilateral
        return [(sigma_1, M), (sigma_1, sigma_n), (m, M), (m, sigma_n)]
    end
end

## -------------------------- band polygons in (C_2, C_3) --------------
"""
    polygon_steady_state(k_1, k_2)

Polygon of the steady-state theorem for the band k_1..k_2, as a LazySets
VPolygon in (C_2, C_3).  Its vertices are the points V_xy at which
g(mu) = det L(mu) has roots x, y and z = S_phi - x - y, for the pairs
(x, y) from `vertex_pairs`.  Returns `nothing` if the polygon is empty or
has zero area.
"""
function polygon_steady_state(k_1, k_2)
    m, sigma_1, sigma_n, M = band_modes(k_1, k_2)
    pairs = vertex_pairs(S_phi, m, sigma_1, sigma_n, M)
    isempty(pairs) && return nothing

    xi2 = phi_C * phi_R + phi_C * phi_rho + phi_R * phi_rho
    xi3 = phi_C * phi_R * phi_rho

    vertices = Vector{Float64}[]
    for (x, y) in pairs
        z   = S_phi - x - y
        C_2 = d_C * d_rho * (xi2 - x * y - (x + y) * z)
        C_3 = f_RR * C_2 - d_C * d_R * d_rho * (x * y * z + xi3)
        push!(vertices, [C_2, C_3])
    end

    polygon = VPolygon(vertices)        # LazySets sorts the vertices itself
    area(polygon) > 0 || return nothing
    return polygon
end

"""
    polygon_hopf(k_1, k_2)

Polygon of the Hopf theorem for the band k_1..k_2, as a LazySets VPolygon
in (C_2, C_3).  Its vertices are the points V̄_xy at which
Q(mu) = a_1 a_2 - a_3 has roots x, y and z = S_theta - x - y, for the
pairs (x, y) from `vertex_pairs`.  Returns `nothing` if the polygon is
empty or has zero area.
"""
function polygon_hopf(k_1, k_2)
    m, sigma_1, sigma_n, M = band_modes(k_1, k_2)
    pairs = vertex_pairs(S_theta, m, sigma_1, sigma_n, M)
    isempty(pairs) && return nothing

    zeta2 = theta_CR * theta_Crho + theta_CR * theta_Rrho + theta_Crho * theta_Rrho
    zeta3 = theta_CR * theta_Crho * theta_Rrho

    p_CR   = d_C + d_R
    p_Crho = d_C + d_rho
    p_Rrho = d_R + d_rho

    vertices = Vector{Float64}[]
    for (x, y) in pairs
        z   = S_theta - x - y
        C_2 = p_CR * p_Rrho * (zeta2 - x * y - (x + y) * z)
        C_3 = p_CR * p_Crho * p_Rrho * (x * y * z + zeta3) - (f_CC + f_rhorho) * C_2
        push!(vertices, [C_2, C_3])
    end

    polygon = VPolygon(vertices)        # LazySets sorts the vertices itself
    area(polygon) > 0 || return nothing
    return polygon
end

## -------------------------- stable zero mode ------------------------
# The zero mode is stable on {g(0) < 0} ∩ {Q(0) > 0}, the "wedge".  Both
# conditions are half-planes in (C_2, C_3).  A LazySets HalfSpace(a, b) is
# the set a . (C_2, C_3) <= b.

# Half-plane g(0) < 0, bounded by the steady-state crossing g(0) = 0:
#     g(0) = f_CC f_RR f_rhorho - f_RR C_2 + C_3 < 0
#         <=>  -f_RR C_2 + C_3 < -f_CC f_RR f_rhorho
const g0_negative = HalfSpace([-f_RR, 1.0], -f_CC * f_RR * f_rhorho)

# Half-plane Q(0) > 0, bounded by the Hopf crossing Q(0) = 0:
#     Q(0) = -(f_CC + f_RR)(f_CC + f_rhorho)(f_RR + f_rhorho) + (f_CC + f_rhorho) C_2 + C_3 > 0
#         <=>  -(f_CC + f_rhorho) C_2 - C_3 < -(f_CC + f_RR)(f_CC + f_rhorho)(f_RR + f_rhorho)
const Q0_positive = HalfSpace([-(f_CC + f_rhorho), -1.0],
                              -(f_CC + f_RR) * (f_CC + f_rhorho) * (f_RR + f_rhorho))

"""
    intersect_zero_mode_stable(polygon)

Intersection of `polygon` with the wedge where the zero mode is stable,
`g0_negative` ∩ `Q0_positive`.  The zero mode also needs
tr L(0) = f_CC + f_RR + f_rhorho < 0, which does not depend on (C_2, C_3).
Returns `nothing` if `polygon` is `nothing` or the intersection is empty
or has zero area.
"""
function intersect_zero_mode_stable(polygon)
    polygon === nothing && return nothing

    region = intersection(intersection(polygon, g0_negative), Q0_positive)
    isempty(region) && return nothing
    area(region) > 0 || return nothing
    return region
end

## -------------------------- no steady-state mode --------------------
"""
    intersect_no_steady_state_mode(region; k_min = κ_min, k_max = κ_max)

Intersection of `region` with {g(-k^2) < 0 for k = k_min..k_max}, the Hopf
theorem's condition g(eta) < 0 checked up to the mode k_max.  For a region
that already has g(0) < 0, this removes the points that lie in a
steady-state polygon of a band within 1..k_max, since the modes with
g > 0 always form a band.

At the mode mu = -k^2,

    g(mu) = (f_CC + mu d_C)(f_RR + mu d_R)(f_rhorho + mu d_rho) - (f_RR + mu d_R) C_2 + C_3,

so g(mu) < 0 is the half-plane  -(f_RR + mu d_R) C_2 + C_3 < -(f_CC + mu d_C)(f_RR + mu d_R)(f_rhorho + mu d_rho).

Modes above k_max are not checked yet, so a point with g > 0 only at
higher modes is kept.

Returns `nothing` if `region` is `nothing` or the intersection is empty or
has zero area.
"""
function intersect_no_steady_state_mode(region; k_min = κ_min, k_max = κ_max)
    region === nothing && return nothing

    for k in k_min:k_max
        mu = -k^2
        g_negative = HalfSpace([-(f_RR + mu * d_R), 1.0],
                               -(f_CC + mu * d_C) * (f_RR + mu * d_R) * (f_rhorho + mu * d_rho))
        region = intersection(region, g_negative)
        isempty(region) && return nothing
    end

    area(region) > 0 || return nothing
    return region
end

## -------------------------- no Hopf mode ----------------------------
"""
    intersect_no_hopf_mode(region; k_min = κ_min, k_max = κ_max)

Intersection of `region` with {Q(-k^2) > 0 for k = k_min..k_max}.  This is
the set P_g'(B) = P_g(B) minus every Hopf polygon P_Q(B') from the
classification theorem, checked up to the mode k_max.  For a region that
already has Q(0) > 0, this removes the points that lie in a Hopf polygon
of a band within 1..k_max, since the modes with Q < 0 always form a band.

At the mode mu = -k^2, write p_ij = d_i + d_j and

    A_ij(mu) = f_ii + f_jj + mu p_ij,

the trace of the (i, j) block of L(mu).  Then

    Q(mu) = -A_CR(mu) A_Crho(mu) A_Rrho(mu) + A_Crho(mu) C_2 + C_3,

so Q(mu) > 0 is the half-plane  -A_Crho(mu) C_2 - C_3 < -A_CR(mu) A_Crho(mu) A_Rrho(mu).
At mu = 0 this is `Q0_positive`.

Modes above k_max are not checked yet, so a point with Q < 0 only at
higher modes is kept.

Returns `nothing` if `region` is `nothing` or the intersection is empty or
has zero area.
"""
function intersect_no_hopf_mode(region; k_min = κ_min, k_max = κ_max)
    region === nothing && return nothing

    for k in k_min:k_max
        mu = -k^2
        A_CR   = f_CC + f_RR     + mu * (d_C + d_R)
        A_Crho = f_CC + f_rhorho + mu * (d_C + d_rho)
        A_Rrho = f_RR + f_rhorho + mu * (d_R + d_rho)
        Q_positive = HalfSpace([-A_Crho, -1.0], -A_CR * A_Crho * A_Rrho)
        region = intersection(region, Q_positive)
        isempty(region) && return nothing
    end

    area(region) > 0 || return nothing
    return region
end

## -------------------------- all bands -------------------------------
"""
    all_bands(k_min = κ_min, k_max = κ_max)

All bands k_1..k_2 with k_min <= k_1 <= k_2 <= k_max, as a vector of
pairs (k_1, k_2).
"""
function all_bands(k_min = κ_min, k_max = κ_max)
    bands = Tuple{Int,Int}[]
    for k_1 in k_min:k_max
        for k_2 in k_1:k_max
            push!(bands, (k_1, k_2))
        end
    end
    return bands
end

## -------------------------- the classification's sets ---------------
# The sets of the classification theorem, each cut to the wedge where the
# zero mode is stable.  A band B is a pair (k_1, k_2).  Each function
# returns a LazySets set, or `nothing` if the set is empty or has zero area.

"""
    steady_state_region(B)

The set P_g'(B): the steady-state polygon of the band B, without the
points that lie in any Hopf polygon.  For every mode in B, L has exactly
one eigenvalue with positive real part, which is real; every other mode
is stable.
"""
function steady_state_region(B)
    polygon = polygon_steady_state(B...)
    return intersect_no_hopf_mode(intersect_zero_mode_stable(polygon))
end

"""
    hopf_region(B)

The set P_Q'(B): the Hopf polygon of the band B, without the points that
lie in any steady-state polygon.  For every mode in B, L has exactly two
eigenvalues with positive real part; every other mode is stable.
"""
function hopf_region(B)
    polygon = polygon_hopf(B...)
    return intersect_no_steady_state_mode(intersect_zero_mode_stable(polygon))
end

"""
    mixed_region(B, B̄)

The set P_gQ(B, B̄): the overlap of the steady-state polygon of the band B
and the Hopf polygon of the band B̄.  For every mode in B, L has exactly
one eigenvalue with positive real part, which is real; for every mode in
B̄ but not in B, exactly two; every other mode is stable.
"""
function mixed_region(B, B̄)
    polygon_g = polygon_steady_state(B...)
    polygon_Q = polygon_hopf(B̄...)
    (polygon_g === nothing || polygon_Q === nothing) && return nothing

    overlap = intersection(polygon_g, polygon_Q)
    isempty(overlap) && return nothing
    return intersect_zero_mode_stable(overlap)
end

"""
    all_regions()

Every non-empty set of the classification theorem, over all bands from
`all_bands`, as a vector of named tuples

    kind      :steady_state for P_g'(B), :hopf for P_Q'(B), :mixed for P_gQ(B, B̄)
    band_g    the band B of the steady-state polygon (`nothing` for :hopf)
    band_Q    the band of the Hopf polygon: B for :hopf, B̄ for :mixed
              (`nothing` for :steady_state)
    region    the LazySets set

By the theorem these sets are mutually disjoint.
"""
function all_regions()
    regions = []
    for B in all_bands() #κ_min and κ_max are specified constants at top of file
        region = steady_state_region(B)
        if region !== nothing
            push!(regions, (kind = :steady_state, band_g = B, band_Q = nothing, region = region))
        end

        region = hopf_region(B)
        if region !== nothing
            push!(regions, (kind = :hopf, band_g = nothing, band_Q = B, region = region))
        end

        for B̄ in all_bands()
            region = mixed_region(B, B̄)
            if region !== nothing
                push!(regions, (kind = :mixed, band_g = B, band_Q = B̄, region = region))
            end
        end
    end
    return regions
end

## -------------------------- sampling the regions --------------------
"""
    sample_region(r, n_points; rng)

`n_points` uniform random points (C_2, C_3) in the region `r.region` of a
labelled region `r` from `all_regions`, as a vector of [C_2, C_3].

Points are drawn uniformly from the region's bounding box and kept if
they lie in the region (LazySets' RejectionSampler).  LazySets' default
sampler is not used: it can switch to a non-uniform method that ignores
the seed.
"""
function sample_region(r, n_points; rng)
    sampler = LazySets.RejectionSampler(r.region)
    return LazySets.sample(r.region, n_points; sampler = sampler, rng = rng)
end

## -------------------------- lift to model rates ---------------------
"""
    lift(C_2, C_3)

Rates `p` for `reaction` such that u* = (1, 1, 1) is a homogeneous steady
state whose Jacobian has the cycle values (C_2, C_3).

We set the off-diagonal entries of the Jacobian according to:

    f_Crho = -1                       (Rho -| Cdc42)
    f_rhoC = -C_2                     (Cdc42 -| Rho, so C_2 = f_Crho f_rhoC)
    f_RC   = 1                        (fixed constant)
    f_rhoR = C_3 / (f_RC f_Crho)      (so C_3 = f_RC f_rhoR f_Crho, i.e. f_rhoR = -C_3)

This matches the signs assumed in cdc42_rac_rho.jl (Cdc42 and Rho inhibit
each other, Rac activates Rho) when C_2 > 0 and C_3 < 0.  For other signs
a Hill amplitude a is negative, which turns that Hill term's inhibition
into activation or the reverse.

We need to choose constants so that the above is true at u = 1.
Looking at the model equation, one choice is the following.

    a_C    = f_Crho / H-'(1; K_C)          a_R    = f_RC   / H+'(1; K_R)
    a_rho1 = f_rhoC / H-'(1; K_rho1)       a_rho2 = f_rhoR / H+'(1; K_rho2)

The diagonal is f_ii = -g_i, from the fixed decay rates.  The basal rates
b_i do not change the Jacobian; they are chosen so that F(1, 1, 1) = 0:

    b_C   = g_C   - a_C H-(1; K_C)
    b_R   = g_R   - a_R H+(1; K_R)
    b_rho = g_rho - a_rho1 H-(1; K_rho1) - a_rho2 H+(1; K_rho2)
"""
function lift(C_2, C_3)
    # Jacobian entries
    f_Crho = -1.0
    f_rhoC = -C_2
    f_rhoR = C_3 / (f_RC * f_Crho)

    # Hill amplitudes from the Hill derivatives at u = 1
    a_C    = f_Crho / dhill_inhibit(1.0, K_C)
    a_R    = f_RC   / dhill_activate(1.0, K_R)
    a_rho1 = f_rhoC / dhill_inhibit(1.0, K_rho1)
    a_rho2 = f_rhoR / dhill_activate(1.0, K_rho2)

    # Basal rates that make u* = (1, 1, 1) a steady state
    b_C   = g_C   - a_C    * hill_inhibit(1.0, K_C)
    b_R   = g_R   - a_R    * hill_activate(1.0, K_R)
    b_rho = g_rho - a_rho1 * hill_inhibit(1.0, K_rho1) - a_rho2 * hill_activate(1.0, K_rho2)

    #In the correct for of p necessary for functions in cdc42_rac_rho.jl 
    return (a_C = a_C, a_R = a_R, a_rho1 = a_rho1, a_rho2 = a_rho2,
            b_C = b_C, b_R = b_R, b_rho = b_rho)
end


## -------------------------- check (1): unstable modes ---------------
"""
    unstable_modes(p, u_star; k_max = 2 * κ_max)

The modes k = 0..k_max that are unstable at the homogeneous steady state
u*: those for which J - k^2 D has an eigenvalue with positive real part.
Here J is the Jacobian of `reaction` at u* and D = diag(d_C, d_R, d_rho).
Modes beyond κ_max are included so that a stray unstable mode would show.

J is computed by automatic differentiation of `reaction` (ForwardDiff),
not from the formulas in `lift`, so this also checks the lift.

Returns a vector with one named tuple per unstable mode:
    k            the mode
    n_positive   the number of eigenvalues with positive real part
    n_complex    how many of those have nonzero imaginary part
    eigenvalues  the eigenvalues with positive real part
"""
function unstable_modes(p, u_star; k_max = 2 * κ_max)
    J = ForwardDiff.jacobian(u -> reaction(u, p), u_star)
    D = Diagonal([d_C, d_R, d_rho])

    unstable = []
    for k in 0:k_max
        positive = filter(λ -> real(λ) > 0, eigvals(J - k^2 * D))
        if !isempty(positive)
            push!(unstable, (k = k, n_positive = length(positive),
                             n_complex = count(λ -> imag(λ) != 0, positive),
                             eigenvalues = positive))
        end
    end
    return unstable
end

"""
    predicted_n_positive(r, k)

The number of eigenvalues of L(-k^2) with positive real part that the
classification theorem predicts at mode k, for a point in the region `r`
from `all_regions`:

    1   if k is in the steady-state band band_g
    2   if k is in the Hopf band band_Q but not in band_g
    0   otherwise (the mode is stable)
"""
function predicted_n_positive(r, k)
    in_band(B) = B !== nothing && B[1] <= k <= B[2]
    in_band(r.band_g) && return 1
    in_band(r.band_Q) && return 2
    return 0
end

"""
    stability_matches(r, modes; k_max = 2 * κ_max)

Whether the unstable modes `modes` (from `unstable_modes`) of a point in
the region `r` match the theorem: for every mode k = 0..k_max, the number
of eigenvalues with positive real part equals `predicted_n_positive`.
"""
function stability_matches(r, modes; k_max = 2 * κ_max)
    observed = Dict(mode.k => mode.n_positive for mode in modes)
    for k in 0:k_max
        get(observed, k, 0) == predicted_n_positive(r, k) || return false
    end
    return true
end

"""
    check_regions(n_points; seed = 1)

For every region from `all_regions`, sample `n_points` points (C_2, C_3),
lift each to model rates, and compare the unstable modes at the
homogeneous steady state with the theorem's prediction for the region.
No PDE is solved.

Returns a vector with one named tuple per point:
    kind       :steady_state, :hopf or :mixed (the region's label)
    band_g     the band of the steady-state polygon, or `nothing`
    band_Q     the band of the Hopf polygon, or `nothing`
    C_2, C_3   the sampled point
    p          the lifted rates
    modes      the unstable modes at u*, from `unstable_modes`
    matches    whether `modes` match the theorem (`stability_matches`)
"""
function check_regions(n_points; seed = 1)
    rng = Xoshiro(seed)
    results = []

    for r in all_regions()
        for (C_2, C_3) in sample_region(r, n_points; rng = rng)
            p       = lift(C_2, C_3)
            u_star  = homogeneous_state(p) #u* ≡ 1 by choice of p
            modes   = unstable_modes(p, u_star)
            matches = stability_matches(r, modes)

            push!(results, (kind = r.kind, band_g = r.band_g, band_Q = r.band_Q,
                            C_2 = C_2, C_3 = C_3,
                            p = p, modes = modes, matches = matches))
        end
    end
    return results
end

"""
    proportion_matching(results)

The proportion of the results from `check_regions` whose stability
matches the theorem.
"""
proportion_matching(results) = count(result -> result.matches, results) / length(results)

## -------------------------- simulating the regions ------------------
"""
    noisy_initial_value(u_star, N; amplitude = 0.01, rng)

Initial value for the PDE: the homogeneous steady state u* on N grid
points, with independent relative noise of size `amplitude` at every
point, u_i(x) = u*_i (1 + amplitude * xi) with xi standard normal.
Returns an N x 3 matrix with columns C, R, rho.
"""
function noisy_initial_value(u_star, N; amplitude = 0.01, rng)
    u_init = zeros(N, 3)
    for species in 1:3
        u_init[:, species] = u_star[species] .* (1 .+ amplitude .* randn(rng, N))
    end
    return u_init
end

"""
    simulate_regions(n_points; N = 257, seed = 1)

For every region from `all_regions`, sample `n_points` points (C_2, C_3),
lift each to model rates, and solve the PDE to steady state from a noisy
initial value around the homogeneous steady state.

Returns a vector with one named tuple per point:
    kind       :steady_state, :hopf or :mixed (the region's label)
    band_g     the band of the steady-state polygon, or `nothing`
    band_Q     the band of the Hopf polygon, or `nothing`
    C_2, C_3   the sampled point
    p          the lifted rates
    u_star     the homogeneous steady state
    u_final    the final N x 3 matrix (columns C, R, rho)
    reached    whether the solver reached a steady state
"""
function simulate_regions(n_points; N = 257, seed = 1)
    rng = Xoshiro(seed)
    results = []

    for r in all_regions()
        for (C_2, C_3) in sample_region(r, n_points; rng = rng)
            p       = lift(C_2, C_3)
            u_star  = homogeneous_state(p) #u* ≡ 1 by choice of p
            u_init  = noisy_initial_value(u_star, N; rng = rng)
            u_final, reached = solve_to_steady_state(u_init, p)

            push!(results, (kind = r.kind, band_g = r.band_g, band_Q = r.band_Q,
                            C_2 = C_2, C_3 = C_3,
                            p = p, u_star = u_star, u_final = u_final, reached = reached))
        end
    end
    return results
end

## -------------------------- check (2): Fourier modes of the final state
"""
    fourier_modes(result; k_max = 2 * κ_max)

Cosine coefficients of the final state of a result from `simulate_regions`,
for the modes k = 0..k_max.  Returns `nothing` if the solver did not reach
a steady state.

With no-flux boundaries on [0, L] = [0, pi], each species expands as

    u(x) = a_0 + sum over k >= 1 of a_k cos(k x),

where

    a_0 = (1 / L) ∫ u(x) dx              (the mean)
    a_k = (2 / L) ∫ u(x) cos(k x) dx     (k >= 1).

The integrals are computed with the trapezoidal rule on the solver's grid
of N equally spaced points from 0 to L.

Returns a vector with one named tuple (k, C, R, rho) per mode, where C, R
and rho are the coefficients a_k of each species.
"""
function fourier_modes(result; k_max = 2 * κ_max)
    result.reached || return nothing

    u = result.u_final
    N = size(u, 1)
    x = range(0, L; length = N)
    h = L / (N - 1)

    # trapezoidal rule: weight h, halved at the two end points
    weights = fill(h, N)
    weights[1]   = h / 2
    weights[end] = h / 2
    integral(v) = sum(weights .* v)

    modes = []
    for k in 0:k_max
        scale = (k == 0) ? 1 / L : 2 / L
        a_k   = [scale * integral(u[:, species] .* cos.(k .* x)) for species in 1:3]
        push!(modes, (k = k, C = a_k[1], R = a_k[2], rho = a_k[3]))
    end
    return modes
end

"""
    dominant_mode(result; k_max = 2 * κ_max, tol = 1e-4, runner_up_ratio = 0.5)

The dominant Fourier mode of the final state of a result from
`simulate_regions`: the mode k = 1..k_max with the largest amplitude

    A_k = sqrt(a_k(C)^2 + a_k(R)^2 + a_k(rho)^2),

combining the three species' coefficients from `fourier_modes`.  The
species are on comparable scales because the lift puts u* at (1, 1, 1).

Returns a named tuple
    k           the dominant mode, or 0 if every A_k is below `tol`
                (the final state is homogeneous)
    runner_up   the mode with the second largest amplitude, if that
                amplitude is at least `runner_up_ratio` times the largest;
                otherwise `nothing`

or `nothing` if the solver did not reach a steady state.
"""
function dominant_mode(result; k_max = 2 * κ_max, tol = 1e-4, runner_up_ratio = 0.5)
    modes = fourier_modes(result; k_max = k_max)
    modes === nothing && return nothing

    # leave out k = 0, which is the mean of each species
    patterned  = filter(mode -> mode.k >= 1, modes)
    amplitudes = [sqrt(mode.C^2 + mode.R^2 + mode.rho^2) for mode in patterned]

    maximum(amplitudes) < tol && return (k = 0, runner_up = nothing)

    # modes sorted from largest to smallest amplitude
    order  = sortperm(amplitudes; rev = true)
    first  = order[1]
    second = order[2]

    k         = patterned[first].k
    # runner_up is the second largest mode if it is at least 0.5 (default) times
    # the largest. Otherwise, it's nothing.
    runner_up = amplitudes[second] >= runner_up_ratio * amplitudes[first] ? patterned[second].k : nothing
    return (k = k, runner_up = runner_up)
end

## -------------------------- plots -----------------------------------
"""
    region_label(r)

Text label for a region `r` from `all_regions`: its band "k_1..k_2" for a
steady-state or Hopf region, and "g k_1..k_2 / Q k̄_1..k̄_2" (the bands of
the steady-state and Hopf polygons) for a mixed region.
"""
function region_label(r)
    band_text(B) = "$(B[1])..$(B[2])"
    r.kind == :steady_state && return band_text(r.band_g)
    r.kind == :hopf         && return band_text(r.band_Q)
    return "g " * band_text(r.band_g) * " / Q " * band_text(r.band_Q)
end

"""
    label_region!(plt, region, label)

Write `label` at the centre of `region` (the average of its vertices,
which lies inside a convex polygon).
"""
function label_region!(plt, region, label)
    vertices = vertices_list(region)
    centre   = sum(vertices) / length(vertices)
    annotate!(plt, centre[1], centre[2], text(label, 8))
end

# Colour and legend entry of each kind of region
const region_colour = Dict(:steady_state => :steelblue, :hopf => :darkorange, :mixed => :seagreen)
const region_legend = Dict(:steady_state => "steady state", :hopf => "Hopf", :mixed => "mixed")

"""
    plot_regions()

All regions from `all_regions` in the (C_2, C_3) plane, coloured by kind
(steady state, Hopf, mixed) and labelled with their bands (`region_label`),
together with the edges of the wedge where the zero mode is stable: the
steady-state crossing g(0) = 0 (solid) and the Hopf crossing Q(0) = 0
(dotted).

With `labels = false` the bands are not written on the regions, which is
useful when zooming in: labels outside the view are still drawn.
"""
function plot_regions(; labels = true)
    plt = plot(xlabel = "C₂", ylabel = "C₃", title = "Regions in cycle space",
               legend = :outerright, size = (900, 500), left_margin = 5Plots.mm)

    # label each kind once in the legend
    labelled = Set{Symbol}()

    for r in all_regions()
        plot!(plt, r.region; color = region_colour[r.kind], alpha = 0.5,
              label = r.kind in labelled ? "" : region_legend[r.kind])
        labels && label_region!(plt, r.region, region_label(r))
        push!(labelled, r.kind)
    end

    # widen the view by 20% on each side, so the wedge edges show beyond the regions
    (x_low, x_high), (y_low, y_high) = xlims(plt), ylims(plt)
    xlims!(plt, x_low - 0.2 * (x_high - x_low), x_high + 0.2 * (x_high - x_low))
    ylims!(plt, y_low - 0.2 * (y_high - y_low), y_high + 0.2 * (y_high - y_low))

    # edges of the wedge where the zero mode is stable
    plot!(plt, Line2D(g0_negative.a, g0_negative.b); color = :black, linealpha = 1,
          linestyle = :solid, label = "g(0) = 0  (steady-state crossing)")
    plot!(plt, Line2D(Q0_positive.a, Q0_positive.b); color = :black, linealpha = 1,
          linestyle = :dot, label = "Q(0) = 0  (Hopf crossing)")
    return plt
end

## -------------------------- real or complex Hopf eigenvalues --------
"""
    plot_hopf_eigenvalues(results)

One panel per Hopf region whose band is a single mode k..k.  Each panel is
the plot of all regions (`plot_regions`) with the sampled points of that
Hopf region, coloured by whether their two unstable eigenvalues at mode k
are a complex pair or two real eigenvalues.

`results` is the output of `check_regions`: one named tuple per sampled
point, with the region's kind and bands, the point (C_2, C_3) and its
unstable `modes` from `unstable_modes`.  The colour is read from
`n_complex` of the mode k.  The output of `simulate_regions` does not work
here, since it has no `modes`.

    results = check_regions(n_points)
    plot_hopf_eigenvalues(results)
"""
function plot_hopf_eigenvalues(results)
    panels = []

    for r in all_regions()
        # only Hopf regions whose band is a single mode
        r.kind == :hopf || continue
        k_1, k_2 = r.band_Q
        k_1 == k_2 || continue
        k = k_1

        plt = plot_regions(labels = false)
        title!(plt, "Hopf region $k..$k")

        # sampled points of this region, split by the eigenvalues at mode k
        complex_C_2, complex_C_3 = Float64[], Float64[]
        real_C_2,    real_C_3    = Float64[], Float64[]
        for result in results
            (result.kind == :hopf && result.band_Q == r.band_Q) || continue
            mode = only(filter(mode -> mode.k == k, result.modes))
            if mode.n_complex == 2
                push!(complex_C_2, result.C_2)
                push!(complex_C_3, result.C_3)
            else
                push!(real_C_2, result.C_2)
                push!(real_C_3, result.C_3)
            end
        end
        scatter!(plt, complex_C_2, complex_C_3; color = :purple, markersize = 3,
                 markerstrokewidth = 0, label = "complex pair")
        scatter!(plt, real_C_2, real_C_3; color = :red, markersize = 3,
                 markerstrokewidth = 0, label = "two real")

        # zoom on the sampled points, with 20% of their spread on each side
        all_C_2 = vcat(complex_C_2, real_C_2)
        all_C_3 = vcat(complex_C_3, real_C_3)
        C_2_low, C_2_high = extrema(all_C_2)
        C_3_low, C_3_high = extrema(all_C_3)
        xlims!(plt, C_2_low - 0.2 * (C_2_high - C_2_low), C_2_high + 0.2 * (C_2_high - C_2_low))
        ylims!(plt, C_3_low - 0.2 * (C_3_high - C_3_low), C_3_high + 0.2 * (C_3_high - C_3_low))

        push!(panels, plt)
    end

    return plot(panels...; layout = (length(panels), 1), size = (900, 500 * length(panels)))
end

## -------------------------- dominant mode of the final state --------
"""
    plot_dominant_modes(results)

The plot of all regions (`plot_regions`) with the sampled points from
`simulate_regions`, coloured by the dominant Fourier mode of their final
state (`dominant_mode`), with one legend entry per mode:

  - filled marker   one mode dominates
  - hollow marker   there is a runner-up mode of similar size
  - grey cross      the solver did not reach a steady state

    results = simulate_regions(n_points)
    plot_dominant_modes(results)
"""
function plot_dominant_modes(results)
    plt = plot_regions()
    title!(plt, "Dominant mode of the final state")

    dominant = [dominant_mode(result) for result in results]

    # points that did not reach a steady state
    not_reached = findall(isnothing, dominant)
    scatter!(plt, [results[i].C_2 for i in not_reached], [results[i].C_3 for i in not_reached];
             marker = :xcross, color = :grey, markersize = 4, label = "no steady state")

    # one colour per dominant mode; mode 0 is the homogeneous state
    reached = filter(!isnothing, dominant)
    for k in sort(unique([d.k for d in reached]))
        colour = (k == 0) ? :black : palette(:tab10)[mod1(k, 10)]
        label  = (k == 0) ? "homogeneous" : "mode $k"

        single = [i for i in eachindex(results) if dominant[i] !== nothing &&
                  dominant[i].k == k && dominant[i].runner_up === nothing]
        mixed  = [i for i in eachindex(results) if dominant[i] !== nothing &&
                  dominant[i].k == k && dominant[i].runner_up !== nothing]

        # filled: one mode dominates
        scatter!(plt, [results[i].C_2 for i in single], [results[i].C_3 for i in single];
                 color = colour, markersize = 4, markerstrokewidth = 0, label = label)
        # hollow: a runner-up mode of similar size
        scatter!(plt, [results[i].C_2 for i in mixed], [results[i].C_3 for i in mixed];
                 markercolor = :transparent, markerstrokecolor = colour,
                 markerstrokewidth = 1.5, markersize = 4, label = "")
    end

    # legend entry explaining the hollow markers
    scatter!(plt, [NaN], [NaN]; markercolor = :transparent, markerstrokecolor = :grey,
             markerstrokewidth = 1.5, label = "hollow: runner-up mode")
    return plt
end




## -------------------------- run -------------------------------------
plot_regions()
plot_hopf_eigenvalues(check_regions(50))
plot_dominant_modes(simulte_regions(50))

