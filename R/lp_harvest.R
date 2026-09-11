# ---------------------------------------------------------------------------
# Section 2.3: the three patterns of exploitation.
#
#   fixed    F_i(t) = F                       Eq. (2.7)
#   BH_P     F_i(t) = c_P   P_i(t)            Eq. (2.8)
#   BH_P/B   F_i(t) = c_P/B P_i(t) / B_i(t)   Eq. (2.9)
#
# All fish enter one mixed-species fishery at a single body mass w_f and
# everything above it is caught at the same rate, so F is size-independent
# above w_f and the balancing is across species, not across sizes.  B_i and
# P_i are measured over that same harvested range, [w_f, w_max,i].
# ---------------------------------------------------------------------------

source("R/lp_model.R")

## --- Measuring B and P over the harvested range ----------------------------
# w_f = 400 g does not fall on a grid point, so the fishery starts at the
# first grid point at or above it; that same point is used for the
# selectivity, for B_i and for P_i, which keeps Y_i = F_i B_i exact.
lp_wf_idx <- function(params, w_f = LP_FISHING$w_f) {
    which(w(params) >= w_f)[1]
}

# B_i(t) = \int_{w_f}^{w_max,i} b_i(w, t) dw                          Eq. (2.3)
lp_biomass <- function(params, n, jlo = other_params(params)$lp_jf) {
    w <- w(params); dw <- dw(params)
    idx <- jlo:length(w)
    drop(n[, idx, drop = FALSE] %*% (w * dw)[idx])
}

# P_i(t) = \int_{w_f}^{w_max,i} [eps(w) g(w,t)/w] b(w,t) dw
#            + [eps g b]_{w_f} - [eps g b]_{w_max,i}                  Eq. (2.4)
#
# Note that the `g` of Eq. (2.4) is the *absolute* rate of mass intake (g/yr),
# not the mass-specific rate `g` of Eq. (A.2) - that is what makes the units of
# both the integrand and the boundary terms come out as g m^-2 yr^-1.  So
# eps*g = dw/dt = mizer's e_growth, b = w phi = w N, and Eq. (2.4) becomes
#
#   P_i = \int e_growth N dw  +  [w e_growth N]_{w_f}  -  [w e_growth N]_{w_max,i}
#
# i.e. production inside the range, plus the biomass flux growing into it at
# w_f, minus the flux growing out of it at the top.  The flux out vanishes
# because all food goes to reproduction at w_max, but it is kept for exactness.
# With the lower bound at egg size the flux in is the rate of egg-mass
# production, w_egg * R_i(t), as the paper notes.
lp_production <- function(params, n, e_growth,
                          jlo = other_params(params)$lp_jf) {
    w <- w(params); dw <- dw(params)
    nw <- length(w)
    idx <- jlo:nw
    interior <- drop((e_growth[, idx, drop = FALSE] * n[, idx, drop = FALSE]) %*%
                         dw[idx])
    flux_in  <- w[jlo] * e_growth[, jlo] * n[, jlo]
    flux_out <- w[nw]  * e_growth[, nw]  * n[, nw]
    interior + flux_in - flux_out
}

# Y_i(t) = F_i(t) B_i(t)                                              Eq. (2.6)
lp_yield <- function(params, n, e_growth, t = 0) {
    lp_F(params, n, e_growth, t) * lp_biomass(params, n)
}

## --- The fishing rules -----------------------------------------------------
# Returns the vector of species-level fishing mortality rates F_i(t).
lp_F <- function(params, n, e_growth, t = 0) {
    op <- other_params(params)
    if (is.null(op$lp_rule) || op$lp_rule == "none")
        return(rep(0, nrow(n)))
    z <- op$lp_zf                                # z'_i, the intensity factor
    F_i <- switch(op$lp_rule,
        fixed = op$lp_const * z,
        BHP   = op$lp_const * z * lp_production(params, n, e_growth),
        BHPB  = op$lp_const * z * lp_production(params, n, e_growth) /
                                  lp_biomass(params, n),
        stop("unknown fishing rule: ", op$lp_rule))
    F_i[!is.finite(F_i)] <- 0
    pmax(F_i, 0)
}

# mizer computes e_growth before f_mort and hands it to the FMort function,
# which is what makes a rule that responds to the current production rate
# possible: project()'s `effort` argument has to be prescribed in advance and
# cannot see the evolving state.
lpFMort <- function(params, n, n_pp, n_other, t = 0, effort,
                    e_growth, pred_mort, ...) {
    outer(lp_F(params, n, e_growth, t), other_params(params)$lp_sel)
}

## --- Setting up a harvesting scenario --------------------------------------
# rule is one of "none", "fixed", "BHP", "BHPB"; const is F, c_P or c_P/B;
# zf is the per-species intensity factor z'_i (all 1 for Figs 3, 4).
lp_set_fishing <- function(params, rule, const = 0,
                           zf = rep(1, nrow(species_params(params))),
                           w_f = LP_FISHING$w_f) {
    jf <- lp_wf_idx(params, w_f)
    other_params(params)$lp_jf    <- jf
    other_params(params)$lp_sel   <- as.numeric(seq_along(w(params)) >= jf)
    other_params(params)$lp_rule  <- rule
    other_params(params)$lp_const <- const
    other_params(params)$lp_zf    <- zf
    setRateFunction(params, "FMort", "lpFMort")
}

## --- Reading B, P, F, Y back off a finished simulation ---------------------
# Recomputes the rates from each saved state, so it works for any rule.
lp_track <- function(sim) {
    p  <- sim@params
    nt <- dim(N(sim))[1]
    tm <- as.numeric(dimnames(N(sim))[[1]])
    out <- lapply(seq_len(nt), function(k) {
        n    <- N(sim)[k, , ]
        n_pp <- NResource(sim)[k, ]
        eg   <- getEGrowth(p, n = n, n_pp = n_pp)
        F_i  <- lp_F(p, n, eg, tm[k])
        B    <- lp_biomass(p, n)
        data.frame(time = tm[k], species = species_params(p)$species,
                   B = B, P = lp_production(p, n, eg), F = F_i, Y = F_i * B,
                   B_total = lp_biomass(p, n, 1),
                   P_total = lp_production(p, n, eg, 1),
                   row.names = NULL)
    })
    do.call(rbind, out)
}

# Convenience: run one scenario from an assembled, relaxed ecosystem.
lp_harvest <- function(params, rule, const = 0,
                       zf = rep(1, nrow(species_params(params))),
                       t_max = LP_FISHING$t_max, dt = LP_NUMERICS$dt,
                       t_save = 1) {
    p <- lp_set_fishing(params, rule, const, zf)
    project(p, t_max = t_max, dt = dt, t_save = t_save, progress_bar = FALSE)
}
