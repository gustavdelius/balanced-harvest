# ---------------------------------------------------------------------------
# Section 3.6 of docs/robustness.md: switching on the larval competition that
# Appendix B nominates as the model's replacement for a stock-recruitment
# relationship, and which at the published parameter values is inert.
#
# The steady state of Eq. (A.11), divided by the carrying capacity a(x) and
# written in terms of the resource level L = n/a and iota = I/a, is
#
#     iota + r L (1 - L) - d L = 0.
#
# At the published values r reaches 40-220/yr over the larval prey range while
# grazing d is under 1/yr, so L sits at 1 whatever iota is: immigration is not
# what protects the plankton, fast regeneration is.  Only once r is brought
# down to the same order as d does iota decide anything.
#
# Total primary production scales as r_0 * a_0, so that product is held fixed
# across the variants.  Otherwise slowing the plankton down would simply
# starve the system, and the comparison would be about how much food there is
# rather than about how contested it is.
# ---------------------------------------------------------------------------

source("R/lp_assembly.R")

lp_plankton_variant <- function(r0, a0, I0, base = LP_PLANKTON)
    modifyList(base, list(r0 = r0, a0 = a0, I0 = I0))

# r_0 * a_0 = 20000 throughout, matching the published pair (10, 2000).
LP_RESOURCE_VARIANTS <- list(
    published = lp_plankton_variant(10,   2e3, 2e3),   # iota = 1
    slow10    = lp_plankton_variant(1,    2e4, 2e3),   # iota = 0.1
    slow100   = lp_plankton_variant(0.1,  2e5, 2e3),   # iota = 0.01
    slow100_noimm = lp_plankton_variant(0.1, 2e5, 0)   # iota = 0
)

# Resource level by size, and the elasticity of growth to fish abundance,
# which is the thing the larval feedback needs to be non-zero.
lp_resource_diagnostics <- function(params, sizes = c(1e-6, 1e-4, 1e-2, 1)) {
    wf <- w_full(params); w <- w(params)
    lvl <- initialNResource(params) / resource_capacity(params)
    lvl[!is.finite(lvl)] <- NA
    n <- initialN(params)

    # re-solve the plankton fixed point for a doubled fish spectrum
    fp <- function(nn) {
        d <- getResourceMort(params, n = nn, n_pp = initialNResource(params))
        r <- resource_rate(params); a <- resource_capacity(params)
        I <- other_params(params)$lp_immigration
        on <- other_params(params)$lp_resource_idx
        out <- rep(0, length(r))
        k <- r[on] / a[on]; b <- r[on] - d[on]
        out[on] <- (b + sqrt(b^2 + 4 * k * I[on])) / (2 * k)
        out
    }
    g1 <- getEGrowth(params)
    g2 <- getEGrowth(params, n = n * 2, n_pp = fp(n * 2))
    elas <- function(tw) {
        k <- which.min(abs(w - tw))
        mean(log(g2[, k] / g1[, k])) / log(2)
    }
    pp <- sum(resource_rate(params) * wf * initialNResource(params) *
                  dw_full(params))
    c(setNames(sapply(sizes, function(s) lvl[which.min(abs(wf - s))]),
               paste0("level_", format(sizes, scientific = TRUE))),
      egg_elasticity = elas(w[1]),
      elasticity_1g  = elas(1),
      primary_production = pp,
      fish_biomass = sum(n %*% (w * dw(params))))
}
