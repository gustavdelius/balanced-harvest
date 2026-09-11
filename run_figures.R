# Apply the three harvesting regimes and draw Figures 2-6.
# Requires data/ecosystems.rds, written by run_assembly.R.
source("R/lp_figures.R")

dir.create("figures", showWarnings = FALSE)
eco <- readRDS("data/ecosystems.rds")

# The runs below take well over an hour, so each expensive step is cached to
# disk and reused if it is already there.  Delete data/cache/ to force a rerun.
dir.create("data/cache", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/cache", paste0(name, ".rds"))
    if (file.exists(f)) { message("  [cached] ", name); return(readRDS(f)) }
    val <- force(expr)
    saveRDS(val, f)
    val
}

total_yield_at <- function(params, rule, const, zf, dt = LP_NUMERICS$dt) {
    sim <- lp_harvest(params, rule, const, zf, dt = dt)
    p <- sim@params
    k <- dim(N(sim))[1]
    n <- N(sim)[k, , ]
    sum(lp_F(p, n, getEGrowth(p, n = n, n_pp = NResource(sim)[k, ])) *
            lp_biomass(p, n))
}

# "For comparability, the constants c_P, c_P/B [...] were calibrated against a
# fixed fishing mortality rate F_i = F = 0.1 /yr, so that all three regimes
# would generate similar ecosystem biomass yields after 50 years."  Their
# ecosystem gave c_P = 1 m^2/g and c_P/B = 0.25; ours is a different
# assemblage, so the same calibration is redone here.  It is done at
# dt = 0.01 to keep it affordable - it only picks a constant - while every
# reported run uses the paper's dt = 0.002.
calibrate <- function(params, target, rule, interval, zf, fallback) {
    f <- function(l) log(total_yield_at(params, rule, exp(l), zf, dt = 0.01)) -
        log(target)
    lo <- log(interval[1]); hi <- log(interval[2])
    flo <- f(lo); fhi <- f(hi)
    if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0) {
        warning(sprintf(
            "could not bracket the yield target for %s over [%g, %g] (yields %.3g to %.3g vs target %.3g); using the paper's value %g",
            rule, interval[1], interval[2], target * exp(flo), target * exp(fhi),
            target, fallback), call. = FALSE)
        return(fallback)
    }
    exp(uniroot(f, lower = lo, upper = hi, f.lower = flo, f.upper = fhi,
                tol = 0.05)$root)
}

## --- Figs 2, 3, 4 ----------------------------------------------------------
p1 <- eco$eco1$params
zf1 <- rep(1, nrow(species_params(p1)))

message("calibrating c_P and c_P/B against F = ", LP_FISHING$F_base)
consts1 <- cached("consts1", {
    target <- total_yield_at(p1, "fixed", LP_FISHING$F_base, zf1, dt = 0.01)
    list(fixed = LP_FISHING$F_base,
         BHP   = calibrate(p1, target, "BHP",  c(0.02, 50), zf1, LP_FISHING$c_P),
         BHPB  = calibrate(p1, target, "BHPB", c(0.01, 3), zf1, LP_FISHING$c_PB))
})
message(sprintf("  c_P = %.3g m2/g   [paper: %g]", consts1$BHP, LP_FISHING$c_P))
message(sprintf("  c_P/B = %.3g      [paper: %g]", consts1$BHPB, LP_FISHING$c_PB))

message("running the three regimes on eco1 at dt = ", LP_NUMERICS$dt)
tr1 <- setNames(lapply(RULES, function(r)
    cached(paste0("tr1_", r),
           lp_track(lp_harvest(p1, r, consts1[[r]], zf1)))), RULES)

# An unfished control over the same 50 years.  The assembled state is only a
# quasi-equilibrium - the paper says as much - so some of the change in the
# relative-biomass panels is drift rather than fishing, and this measures it.
ctrl1 <- cached("ctrl1", lp_track(lp_harvest(p1, "none", 0, zf1)))

f2 <- lp_fig2(p1)
ggsave("figures/fig2.pdf", f2$plot, width = 7, height = 3)
ggsave("figures/fig3.pdf", lp_fig3(tr1), width = 7.5, height = 7)

# The paper traces one rare species through Fig. 4; take the rarest fished one.
st <- subset(tr1$fixed, time == 0 & B > 0)
rare <- st$species[which.min(st$B)]
f4 <- lp_fig4(tr1, p1, E = consts1$BHPB, track_species = rare)
ggsave("figures/fig4.pdf", f4$plot, width = 7.5, height = 3)

## --- Figs 5, 6 -------------------------------------------------------------
# The paper doubled the baseline and randomised the intensity, keeping the
# constants at exactly twice their Figs 3, 4 values.
consts6 <- lapply(consts1, function(x) 2 * x)
set.seed(99)
zfs <- lapply(eco$ecor, function(e)
    runif(nrow(species_params(e$params)), LP_FISHING$z_range[1],
          LP_FISHING$z_range[2]))

ggsave("figures/fig5.pdf", lp_fig5(lapply(eco$ecor, `[[`, "params")),
       width = 7.5, height = 3)

message("running the three regimes on the three randomised ecosystems")
tr6 <- lapply(seq_along(eco$ecor), function(k)
    setNames(lapply(RULES, function(r)
        cached(sprintf("tr6_%d_%s", k, r),
               lp_track(lp_harvest(eco$ecor[[k]]$params, r, consts6[[r]],
                                   zfs[[k]])))), RULES))
ggsave("figures/fig6.pdf", lp_fig6(tr6), width = 7.5, height = 7)

# Unfished controls for the other three ecosystems too, so that every
# comparison below can be made against no fishing rather than against year 0.
ctrl6 <- lapply(seq_along(eco$ecor), function(k)
    cached(sprintf("ctrl6_%d", k),
           lp_track(lp_harvest(eco$ecor[[k]]$params, "none", 0, zfs[[k]]))))

## --- Numbers the paper quotes ----------------------------------------------
pp_rate <- sum(resource_rate(p1) * w_full(p1) * initialNResource(p1) *
                   dw_full(p1))
# `rel` is biomass at year 50 over biomass at year 0, as the paper reports it.
# `relc` is biomass at year 50 over the unfished control at year 50, which
# separates the effect of fishing from the drift of a quasi-equilibrium.
summarise <- function(tracks, control, label, rules = RULES) {
    cend <- subset(control, time == LP_FISHING$t_max)
    do.call(rbind, lapply(rules, function(r) {
        d <- tracks[[r]]
        e <- subset(d, time == LP_FISHING$t_max)
        s <- subset(d, time == 0)
        rel  <- e$B_total[match(s$species, e$species)] / s$B_total
        relc <- e$B_total[match(cend$species, e$species)] / cend$B_total
        data.frame(ecosystem = label, rule = r,
                   yield = sum(e$Y),
                   fogarty_permille = 1000 * sum(e$Y) / pp_rate,
                   min_rel = min(rel), min_rel_ctrl = min(relc),
                   n_below_10pct = sum(rel < 0.1),
                   n_below_10pct_ctrl = sum(relc < 0.1),
                   n_below_1pct_ctrl  = sum(relc < 0.01))
    }))
}
smry <- rbind(summarise(list(none = ctrl1), ctrl1, "eco1 (unfished control)",
                        "none"),
              summarise(tr1, ctrl1, "eco1"),
              do.call(rbind, lapply(seq_along(tr6), function(k)
                  summarise(tr6[[k]], ctrl6[[k]], paste0("eco_r", k)))))
print(smry, row.names = FALSE, digits = 3)
message("\nB ~ P^alpha over the harvested range: alpha = ", signif(f4$alpha, 4),
        "   [paper: 1.004]")

saveRDS(list(tr1 = tr1, ctrl1 = ctrl1, ctrl6 = ctrl6, tr6 = tr6, zfs = zfs, consts1 = consts1,
             consts6 = consts6, alpha = f4$alpha, fig2 = f2$data,
             summary = smry, pp_rate = pp_rate, rare = rare),
        "data/results.rds")
message("wrote data/results.rds and figures/fig{2,3,4,5,6}.pdf")
