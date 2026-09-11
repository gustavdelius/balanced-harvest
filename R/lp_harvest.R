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

# The three rules of the paper, in the order its figures use.
RULES <- c("fixed", "BHP", "BHPB")

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
# Every rule has the form F_i(t) = c * z'_i * g_i(t), differing only in the
# per-species allocation g_i.  The three of the paper are `fixed`, `BHP` and
# `BHPB`; the rest are controls for docs/robustness.md.
#
#   fixed   g_i = 1                        Eq. (2.7)
#   BHP     g_i = P_i(t)                   Eq. (2.8)
#   BHPB    g_i = P_i(t)/B_i(t)            Eq. (2.9)
#   BHB     g_i = B_i(t)                   does production do any work, or is
#                                          proportionality to biomass enough?
#   power   g_i = B_i(t)^theta             the family the three rules sit in:
#                                          theta = 0 is fixed, 1 is BHB,
#                                          BHPB sits at theta slightly negative
#   frozen  g_i = a fixed vector           separates the allocation across
#                                          species from the feedback over time
# P_i and B_i are measured from `lp_jm`, which the paper sets equal to the
# harvested range (`lp_jf`).  Keeping them separate lets the measurement range
# be varied without changing what is caught (robustness.md item 11).
lp_alloc <- function(params, n, e_growth) {
    op <- other_params(params)
    jm <- if (is.null(op$lp_jm)) op$lp_jf else op$lp_jm
    switch(op$lp_rule,
        none   = rep(0, nrow(n)),
        fixed  = rep(1, nrow(n)),
        BHP    = lp_production(params, n, e_growth, jm),
        BHPB   = lp_production(params, n, e_growth, jm) /
                     lp_biomass(params, n, jm),
        BHB    = lp_biomass(params, n, jm),
        power  = lp_biomass(params, n, jm)^op$lp_theta,
        frozen = op$lp_alloc,
        stop("unknown fishing rule: ", op$lp_rule))
}

# Returns the vector of species-level fishing mortality rates F_i(t).
lp_F <- function(params, n, e_growth, t = 0) {
    op <- other_params(params)
    if (is.null(op$lp_rule) || op$lp_rule == "none")
        return(rep(0, nrow(n)))
    F_i <- op$lp_const * op$lp_zf * lp_alloc(params, n, e_growth)
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
                           w_f = LP_FISHING$w_f,
                           theta = NULL, alloc = NULL, w_measure = NULL) {
    jf <- lp_wf_idx(params, w_f)
    other_params(params)$lp_jf    <- jf
    other_params(params)$lp_jm    <- if (is.null(w_measure)) jf
                                     else lp_wf_idx(params, w_measure)
    other_params(params)$lp_sel   <- as.numeric(seq_along(w(params)) >= jf)
    other_params(params)$lp_rule  <- rule
    other_params(params)$lp_const <- const
    other_params(params)$lp_zf    <- zf
    other_params(params)$lp_theta <- theta
    other_params(params)$lp_alloc <- alloc
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
                       t_save = 1, theta = NULL, alloc = NULL,
                       w_f = LP_FISHING$w_f, w_measure = NULL) {
    p <- lp_set_fishing(params, rule, const, zf, w_f = w_f,
                        theta = theta, alloc = alloc, w_measure = w_measure)
    project(p, t_max = t_max, dt = dt, t_save = t_save, progress_bar = FALSE)
}
