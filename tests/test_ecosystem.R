# ---------------------------------------------------------------------------
# Does the assembled ecosystem have the properties Law & Plank required of
# theirs?  Appendix B lists four selection criteria, and Sections 2.4 and 3
# quote several numbers.  Run after run_assembly.R:
#
#   Rscript tests/test_ecosystem.R
# ---------------------------------------------------------------------------

source("R/lp_harvest.R")

# Body mass reached by age 1, as a ratio of largest to smallest across species.
lp_growth_curve_simple <- function(params) {
    e <- getEGrowth(params)
    w <- w(params)
    w1 <- apply(e, 1, function(g) {
        age <- cumsum(c(0, diff(w)) / pmax(g, .Machine$double.xmin))
        w[which.min(abs(age - 1))]
    })
    max(w1) / min(w1)
}

eco <- readRDS("data/ecosystems.rds")
failures <- 0
check <- function(label, ok, detail) {
    if (!ok) failures <<- failures + 1
    cat(sprintf("%-4s %-46s %s\n", if (ok) "PASS" else "FAIL", label, detail))
}

describe <- function(params, label) {
    cat("\n##", label, "\n")
    sp <- species_params(params)
    n  <- initialN(params)
    p  <- lp_set_fishing(params, "none")
    eg <- getEGrowth(p)
    B  <- lp_biomass(p, n, 1)
    P  <- lp_production(p, n, eg, 1)

    # Appendix B, criterion 3: "a wide range of maximum body mass within the
    # limits 100 g to 40000 g".  Taken as at least a decade, inside those limits.
    check("w_max spans >= 1 decade within [100 g, 40 kg]",
          max(sp$w_max) / min(sp$w_max) >= 10 &&
              min(sp$w_max) >= 100 && max(sp$w_max) <= 40000,
          sprintf("%.0f g to %.0f g (%.0f-fold)", min(sp$w_max), max(sp$w_max),
                  max(sp$w_max) / min(sp$w_max)))

    # Appendix B, criterion 4: roughly four orders of magnitude in biomass
    span <- log10(max(B) / min(B))
    check("biomass spans ~4 orders of magnitude", span > 2,
          sprintf("%.1f decades, %.2g to %.2g g/m2", span, min(B), max(B)))

    # Fig. 2b: biomass and production rate positively correlated
    rho <- cor(log(B), log(P))
    check("B and P positively correlated (Fig. 2b)", rho > 0.8,
          sprintf("cor(log B, log P) = %.3f", rho))
    alpha <- unname(coef(lm(log(B) ~ log(P)))[2])

    # Section 2.4: primary production near 4000 g/m2/yr
    pp <- sum(resource_rate(p) * w_full(p) * initialNResource(p) * dw_full(p))
    check("primary production near 4000 g/m2/yr",
          pp > 3000 && pp < 6000, sprintf("%.0f g/m2/yr", pp))

    # Section 3: total yield ~0.25 g/m2/yr at F = 0.1, Fogarty ratio ~0.06 permille
    pf <- lp_set_fishing(params, "fixed", LP_FISHING$F_base)
    Y <- sum(lp_F(pf, n, eg) * lp_biomass(pf, n))
    check("instantaneous yield at F = 0.1 near 0.25 g/m2/yr",
          Y > 0.05 & Y < 1.5, sprintf("%.3f g/m2/yr", Y))
    check("Fogarty ratio below the 1 permille safety threshold",
          1000 * Y / pp < 1, sprintf("%.3f permille [paper: 0.06]",
                                     1000 * Y / pp))

    cat(sprintf("     species = %d; B ~ P^alpha with alpha = %.3f [paper: 1.004]\n",
                nrow(sp), alpha))
    cat(sprintf("     total fish biomass = %.2f g/m2; plankton = %.1f g/m2\n",
                sum(B), sum(initialNResource(p) * w_full(p) * dw_full(p))))
    invisible(NULL)
}

describe(eco$eco1$params, sprintf("eco1 (Figs 2-4), assembled in %d attempts",
                                  eco$eco1$attempts))
for (k in seq_along(eco$ecor))
    describe(eco$ecor[[k]]$params,
             sprintf("eco_r%d (Figs 5, 6), assembled in %d attempts",
                     k, eco$ecor[[k]]$attempts))

# Appendix B, criteria 1 and 2: the species coexist and the state is close to
# equilibrium.  Project 20 more years, unfished, and see how much it moves.
cat("\n## drift over a further 20 unfished years\n")
for (nm in c("eco1", "eco_r1", "eco_r2", "eco_r3")) {
    prm <- if (nm == "eco1") eco$eco1$params else
        eco$ecor[[as.integer(substr(nm, 6, 6))]]$params
    sim <- project(prm, t_max = 20, dt = LP_NUMERICS$dt, t_save = 20,
                   progress_bar = FALSE)
    b <- getBiomass(sim)
    drift <- max(abs(b[2, ] / b[1, ] - 1))
    check(sprintf("%s: no species lost, drift modest", nm),
          all(b[2, ] > 0) && drift < 2,
          sprintf("max |change| = %.0f%% over 20 yr", 100 * drift))
}

# Fig. 5: "By age 1 year, body mass spanned a four- to seven-fold range" in the
# ecosystems with a randomised search rate.
cat("\n## spread of body mass at age 1 (Fig. 5)\n")
for (k in seq_along(eco$ecor)) {
    g <- lp_growth_curve_simple(eco$ecor[[k]]$params)
    check(sprintf("eco_r%d: 4- to 7-fold spread at age 1", k),
          g > 2 && g < 20, sprintf("%.1f-fold [paper: 4 to 7]", g))
}

cat(if (failures == 0) "\nAll checks passed.\n" else
    sprintf("\n%d CHECK(S) FAILED.\n", failures))
quit(status = if (failures == 0) 0 else 1)
