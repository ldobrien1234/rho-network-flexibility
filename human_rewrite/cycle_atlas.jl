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

# --------------------------- fixed parameters ------------------------
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
const κ_max = 4

# Off-diagonal Jacobian entry Cdc42 -> Rac, held fixed in the lift
const f_RC = 1.0

# --------------------------- vertex pairs ----------------------------
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

# --------------------------- band polygons in (C_2, C_3) --------------
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

# --------------------------- stable zero mode ------------------------
"""
    intersect_zero_mode_stable(polygon)

Intersection of `polygon` with the region where the zero mode is stable,
{g(0) < 0} ∩ {Q(0) > 0}.  Both conditions are half-planes in (C_2, C_3):

    g(0) = f_CC f_RR f_rhorho - f_RR C_2 + C_3 < 0
        <=>  -f_RR C_2 + C_3 < -f_CC f_RR f_rhorho

    Q(0) = -(f_CC + f_RR)(f_CC + f_rhorho)(f_RR + f_rhorho) + (f_CC + f_rhorho) C_2 + C_3 > 0
        <=>  -(f_CC + f_rhorho) C_2 - C_3 < -(f_CC + f_RR)(f_CC + f_rhorho)(f_RR + f_rhorho)

A LazySets HalfSpace(a, b) is the set a . (C_2, C_3) <= b.  The zero mode
also needs tr L(0) = f_CC + f_RR + f_rhorho < 0, which does not depend on
(C_2, C_3).  Returns `nothing` if `polygon` is `nothing` or the
intersection is empty or has zero area.
"""
function intersect_zero_mode_stable(polygon)
    polygon === nothing && return nothing

    g0_negative = HalfSpace([-f_RR, 1.0], -f_CC * f_RR * f_rhorho)
    Q0_positive = HalfSpace([-(f_CC + f_rhorho), -1.0],
                            -(f_CC + f_RR) * (f_CC + f_rhorho) * (f_RR + f_rhorho))

    region = intersection(intersection(polygon, g0_negative), Q0_positive)
    isempty(region) && return nothing
    area(region) > 0 || return nothing
    return region
end

# --------------------------- no steady-state mode --------------------
"""
    intersect_no_steady_state_mode(region; k_max = κ_max)

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
function intersect_no_steady_state_mode(region; k_max = κ_max)
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



# --------------------------- all bands -------------------------------
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

# --------------------------- sampling the band polygons --------------

"""
    sample_band(k_1, k_2, n_points; seed = 1)

Uniform random points (C_2, C_3) in the band k_1..k_2's polygon for the
steady-state theorem (A) and for the Hopf theorem (B), each cut to the
region where the zero mode is stable (`intersect_zero_mode_stable`).  The
Hopf region is also cut to g < 0 at every mode
(`intersect_no_steady_state_mode`), which removes the points that lie in
a steady-state polygon, where the steady-state behaviour dominates.
Returns a named tuple `(steady_state, hopf)` of vectors of [C_2, C_3]; a
vector is empty when its cut region is.

Points are drawn uniformly from the polygon's bounding box and kept if
they lie in the polygon (LazySets' RejectionSampler).  LazySets' default
sampler is not used: it can switch to a non-uniform method that ignores
the seed.
"""
function sample_band(k_1, k_2, n_points; seed = 1)
    rng = Xoshiro(seed)

    region_A = intersect_zero_mode_stable(polygon_steady_state(k_1, k_2))
    region_B = intersect_no_steady_state_mode(intersect_zero_mode_stable(polygon_hopf(k_1, k_2)))

    points_A = region_A === nothing ? Vector{Float64}[] :
        LazySets.sample(region_A, n_points; sampler = LazySets.RejectionSampler(region_A), rng = rng)
    points_B = region_B === nothing ? Vector{Float64}[] :
        LazySets.sample(region_B, n_points; sampler = LazySets.RejectionSampler(region_B), rng = rng)

    return (steady_state = points_A, hopf = points_B)
end

# --------------------------- lift to model rates ---------------------
"""
    lift(C_2, C_3)

Rates `p` for `reaction` such that u* = (1, 1, 1) is a homogeneous steady
state whose Jacobian has the cycle values (C_2, C_3).

We set the off-diagonal entries of the Jacobian according to:

    f_Crho = f_rhoC = -sqrt(C_2)      (mutual inhibition, so C_2 = f_Crho f_rhoC)
    f_RC   = 1                        (fixed constant)
    f_rhoR = C_3 / (f_RC f_Crho)      (so C_3 = f_RC f_rhoR f_Crho)

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
    f_Crho = -sqrt(C_2)
    f_rhoC = -sqrt(C_2)
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


# --------------------------- simulating the bands --------------------
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
    simulate_bands(n_points; N = 257, seed = 1)

For every band, sample `n_points` points (C_2, C_3) from its steady-state
and Hopf regions, lift each to model rates, and solve the PDE to steady
state from a noisy initial value around the homogeneous steady state.

Returns a vector with one named tuple per point:
    band       (k_1, k_2)
    theorem    :steady_state or :hopf
    C_2, C_3   the sampled point
    p          the lifted rates
    u_star     the homogeneous steady state
    u_final    the final N x 3 matrix (columns C, R, rho)
    reached    whether the solver reached a steady state
"""
function simulate_bands(n_points; N = 257, seed = 1)
    rng = Xoshiro(seed)
    results = []

    for (k_1, k_2) in all_bands()
        points = sample_band(k_1, k_2, n_points; seed = seed)

        for theorem in (:steady_state, :hopf)
            for (C_2, C_3) in points[theorem]
                p       = lift(C_2, C_3)
                u_star  = homogeneous_state(p) #u*≡1 by choice of p
                u_init  = noisy_initial_value(u_star, N; rng = rng)
                u_final, reached = solve_to_steady_state(u_init, p)

                push!(results, (band = (k_1, k_2), theorem = theorem, C_2 = C_2, C_3 = C_3,
                                p = p, u_star = u_star, u_final = u_final, reached = reached))
            end
        end
    end
    return results
end
