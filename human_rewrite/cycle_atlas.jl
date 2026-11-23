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

# --------------------------- fixed parameters ------------------------
# Diagonal of the Jacobian.  No species regulates itself, so f_ii comes
# only from linear decay: f_ii = -g_i.
const f_CC     = -g_C
const f_RR     = -g_R
const f_rhorho = -g_rho

# phi_i = f_ii / d_i and S_phi = -(phi_C + phi_R + phi_rho)
const phi_C   = f_CC     / d_C
const phi_R   = f_RR     / d_R
const phi_rho = f_rhorho / d_rho
const S_phi   = -(phi_C + phi_R + phi_rho)

# Bands k_1..k_2 are checked for 1 <= k_1 <= k_2 <= k_max
const k_max = 4

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

The pairs (x, y) of the vertices V_xy of the region for a band, read off
from the table by the value of S = S_phi.  Each row of the table is one
branch below.  The thresholds are in increasing order because
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
