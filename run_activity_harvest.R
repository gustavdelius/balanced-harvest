# The three harvesting rules on an ecosystem with a real fast-slow continuum.
#
# In the paper's parameterisation P/B is nearly the same for every species, so
# BH_P/B - a constant exploitation ratio - is barely distinguishable from a
# constant F.  That is the argument of Fig. 3's middle row.  The neutral-
# activity ecosystems built by run_activity_assembly.R have P/B spread over
# 5-8x instead of 1.7-2x, so the argument no longer applies to them and
# BH_P/B should fish fast species several times harder than slow ones.
#
#   Rscript run_activity_harvest.R <seed>
#
source("R/lp_harvest.R")

seed <- as.integer(commandArgs(trailingOnly = TRUE)[1])
p <- readRDS(sprintf("data/activity/neutral_%d.rds", seed))$params
zf <- rep(1, nrow(species_params(p)))

# The fast end runs ~2.9x the paper's rates, so the reported runs use
# dt = 0.001 rather than the paper's 0.002; calibration, which only picks a
# constant, uses dt = 0.004.  See the numerics section of docs/life-history.md.
DT_RUN <- 0.001
DT_CAL <- 0.004

total_yield_at <- function(rule, const, dt) {
    sim <- lp_harvest(p, rule, const, zf, dt = dt)
    q <- sim@params; k <- dim(N(sim))[1]; n <- N(sim)[k, , ]
    sum(lp_F(q, n, getEGrowth(q, n = n, n_pp = NResource(sim)[k, ])) *
            lp_biomass(q, n))
}

calibrate <- function(target, rule, interval, fallback) {
    f <- function(l) log(total_yield_at(rule, exp(l), DT_CAL)) - log(target)
    lo <- log(interval[1]); hi <- log(interval[2])
    flo <- f(lo); fhi <- f(hi)
    if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
        warning(sprintf("could not bracket %s over [%g, %g]; using %g",
                        rule, interval[1], interval[2], fallback), call. = FALSE)
        return(fallback)
    }
    exp(uniroot(f, lower = lo, upper = hi, f.lower = flo, f.upper = fhi,
                tol = 0.05)$root)
}

message("seed ", seed, ": calibrating against F = ", LP_FISHING$F_base)
target <- total_yield_at("fixed", LP_FISHING$F_base, DT_CAL)
consts <- list(fixed = LP_FISHING$F_base,
               BHP   = calibrate(target, "BHP",  c(0.02, 50), LP_FISHING$c_P),
               BHPB  = calibrate(target, "BHPB", c(0.01, 3), LP_FISHING$c_PB))
message(sprintf("  c_P = %.3g   c_P/B = %.3g", consts$BHP, consts$BHPB))

message("  running the three regimes and an unfished control at dt = ", DT_RUN)
tr <- setNames(lapply(RULES, function(r)
    lp_track(lp_harvest(p, r, consts[[r]], zf, dt = DT_RUN))), RULES)
ctrl <- lp_track(lp_harvest(p, "none", 0, zf, dt = DT_RUN))

dir.create("data/activity_harvest", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(seed = seed, params = p, consts = consts, tracks = tr,
             control = ctrl, z = species_params(p)$z),
        sprintf("data/activity_harvest/harvest_%d.rds", seed))
message("  wrote data/activity_harvest/harvest_", seed, ".rds")
