# ---------------------------------------------------------------------------
# Balanced harvesting in the sense of Law & Plank (2023), "Fishing for
# biodiversity by balanced harvesting", Fish and Fisheries, doi:10.1111/faf.12705
#
#   fixed    F_i(t) = F_i                      (2.7)
#   BH_P     F_i(t) = c_P     P_i(t)           (2.8)
#   BH_P/B   F_i(t) = c_P/B   P_i(t) / B_i(t)  (2.9)
#
# In all three, fishing mortality is independent of body size above a single
# shared entry mass w_f: all fish enter one mixed-species fishery at w_f and
# everything above it is caught at the same rate. The balancing is across
# species, not across sizes.
# ---------------------------------------------------------------------------

library(mizer)

W_F    <- 400   # g, common knife-edge entry mass (the value used in the paper)
T_MAX  <- 50    # years of fishing
F_FIXED <- 0.2  # baseline fishing mortality for the fixed-F regime

## Model -------------------------------------------------------------------
# One gear per species, so that each species' fishing mortality can be set
# independently, with a knife edge at the shared entry mass W_F.
setup <- function(params = NS_params, w_f = W_F) {
    sp <- species_params(params)$species
    gear_params(params) <- data.frame(
        gear = sp, species = sp, sel_func = "knife_edge",
        knife_edge_size = w_f, catchability = 1)
    initial_effort(params) <- setNames(rep(0, length(sp)), sp)
    params
}

## Somatic production rate over the harvested size range --------------------
# P_i = \int_{w_f}^{w_max} g_i(w) n_i(w) dw   (interior term of their Eq. 2.4)
production <- function(params, n, e_growth) {
    drop(sweep(e_growth * n, 2, params@dw, "*") %*% params@selectivity[1, 1, ])
}
biomass_fished <- function(params, n) {
    drop(sweep(n, 2, params@w * params@dw, "*") %*% params@selectivity[1, 1, ])
}

## The fishing rules as a custom FMort --------------------------------------
# mizer computes e_growth before f_mort and passes it in, so the production
# rate is available at the moment the rule needs it.
balancedFMort <- function(params, n, n_pp, n_other, t, effort,
                          e_growth, pred_mort, ...) {
    op  <- other_params(params)
    sel <- params@selectivity[1, 1, ]
    F_i <- switch(op$rule,
        none   = op$z * 0,
        fixed  = op$const * op$z,
        BHP    = op$const * op$z * production(params, n, e_growth),
        BHPB   = op$const * op$z * production(params, n, e_growth) /
                                   biomass_fished(params, n))
    F_i[!is.finite(F_i)] <- 0
    outer(F_i, sel)
}

setRule <- function(params, rule, const, z) {
    other_params(params)$rule  <- rule
    other_params(params)$const <- const
    other_params(params)$z     <- z
    setRateFunction(params, "FMort", "balancedFMort")
}

## Recover F_i(t) from a saved simulation -----------------------------------
FOverTime <- function(sim) {
    p <- sim@params
    t(vapply(seq_len(dim(N(sim))[1]), function(k) {
        n    <- N(sim)[k, , ]
        n_pp <- NResource(sim)[k, ]
        balancedFMort(p, n, n_pp, list(), 0, 0,
                      getEGrowth(p, n = n, n_pp = n_pp), NULL)[, length(p@w)]
    }, numeric(nrow(species_params(p)))))
}
