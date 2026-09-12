# BH_P against BH_B on the fast-slow ecosystems.
#
# robustness-results.md finds that F proportional to biomass does everything
# F proportional to production does.  That finding is a corollary of P/B being
# nearly the same for every species: if P_i = k B_i with k common to all, the
# two allocations are the same up to the calibration constant.  Once P/B varies
# - as it does by 5-8x in the neutral ecosystems of life-history.md - they come
# apart, and by a known amount:
#
#   F_i(BH_P) ~ P_i = (P_i/B_i) B_i ~ z_i B_i      versus     F_i(BH_B) ~ B_i
#
# so BH_P is BH_B tilted by the activity factor.  Which is better depends on
# whether that tilt carries information about which species are at risk.  It
# should not: rarity is what matters and z is orthogonal to it.
#
#   Rscript run_activity_bhb.R <seed>
#
source("R/lp_harvest.R")

seed <- as.integer(commandArgs(trailingOnly = TRUE)[1])
f_out <- sprintf("data/activity_harvest/harvest_%d.rds", seed)
h <- readRDS(f_out)
p <- h$params
zf <- rep(1, nrow(species_params(p)))
DT_RUN <- 0.001
DT_CAL <- 0.004

total_yield_at <- function(rule, const, dt) {
    sim <- lp_harvest(p, rule, const, zf, dt = dt)
    q <- sim@params; k <- dim(N(sim))[1]; n <- N(sim)[k, , ]
    sum(lp_F(q, n, getEGrowth(q, n = n, n_pp = NResource(sim)[k, ])) *
            lp_biomass(q, n))
}

# Recomputed exactly as run_activity_harvest.R did, so the four rules are
# calibrated against one and the same target.
message("seed ", seed, ": recomputing the yield target")
target <- total_yield_at("fixed", LP_FISHING$F_base, DT_CAL)

f <- function(l) log(total_yield_at("BHB", exp(l), DT_CAL)) - log(target)
lo <- log(1e-3); hi <- log(1e3)
flo <- f(lo); fhi <- f(hi)
if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0)
    stop("could not bracket the yield target for BHB")
c_B <- exp(uniroot(f, lower = lo, upper = hi, f.lower = flo, f.upper = fhi,
                   tol = 0.05)$root)
message(sprintf("  c_B = %.3g   (target yield %.4g)", c_B, target))

message("  running BH_B at dt = ", DT_RUN)
h$tracks$BHB <- lp_track(lp_harvest(p, "BHB", c_B, zf, dt = DT_RUN))
h$consts$BHB <- c_B
saveRDS(h, f_out)
message("  updated ", f_out)
