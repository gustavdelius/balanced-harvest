# Phases 1 and 2 of docs/robustness.md.
#
#   Phase 1  the yield-biodiversity frontier of the three published rules, on
#            all four ecosystems, which removes the "equal yield after 50
#            years" calibration choice from the comparison entirely; plus the
#            alternative calibration criteria of item 3.
#   Phase 2  controls that test whether the mechanism is the one claimed:
#            F proportional to biomass rather than production, the frozen
#            allocations that separate allocation from feedback, and the
#            theta family F ~ B^theta that all three rules sit inside.
#
# Sweeps run at dt = 0.01 (checked against the paper's dt = 0.002 in
# tests/test_convergence.R); controls are run at the same dt so that
# numerical differences cancel.
source("R/lp_experiments.R")
library(parallel)
NCORES <- as.integer(Sys.getenv("LP_CORES", "2"))

dir.create("data/robustness", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/robustness", paste0(name, ".rds"))
    if (file.exists(f)) { message("  [cached] ", name); return(readRDS(f)) }
    val <- force(expr); saveRDS(val, f); val
}

eco <- readRDS("data/ecosystems.rds")
set.seed(99)
zfs <- lapply(eco$ecor, function(e)
    runif(nrow(species_params(e$params)), LP_FISHING$z_range[1],
          LP_FISHING$z_range[2]))
ECOS <- list(eco1   = list(p = eco$eco1$params, zf = NULL, F_ref = 0.1),
             eco_r1 = list(p = eco$ecor[[1]]$params, zf = zfs[[1]], F_ref = 0.2),
             eco_r2 = list(p = eco$ecor[[2]]$params, zf = zfs[[2]], F_ref = 0.2),
             eco_r3 = list(p = eco$ecor[[3]]$params, zf = zfs[[3]], F_ref = 0.2))

controls <- lapply(names(ECOS), function(k)
    cached(paste0("ctrl_", k), lp_control(ECOS[[k]]$p, ECOS[[k]]$zf)))
names(controls) <- names(ECOS)

## --- Phase 1: frontiers of the three published rules -----------------------
message("Phase 1: frontiers of fixed, BH_P and BH_P/B")
jobs1 <- expand.grid(eco = names(ECOS), rule = RULES, stringsAsFactors = FALSE)
front1 <- do.call(rbind, mclapply(seq_len(nrow(jobs1)), function(i) {
    k <- jobs1$eco[i]; r <- jobs1$rule[i]; e <- ECOS[[k]]
    cbind(ecosystem = k,
          cached(sprintf("front_%s_%s", k, r),
                 lp_frontier(e$p, r, controls[[k]], e$zf, F_ref = e$F_ref)))
}, mc.cores = NCORES))

## --- Phase 1 item 3: does the calibration criterion change c_P? ------------
message("Phase 1: alternative calibration criteria")
alt_calib <- cached("alt_calib", {
    e <- ECOS$eco1
    # reference yields from the fixed-F run at multiplier 1
    f <- subset(front1, ecosystem == "eco1" & rule == "fixed")
    f1 <- subset(f, mult == 1)
    targets <- c(terminal = f1$yield, cumulative = f1$yield_cum)
    do.call(rbind, lapply(RULES[-1], function(r) {
        fr <- subset(front1, ecosystem == "eco1" & rule == r)
        data.frame(rule = r,
                   c_terminal  = lp_at_yield(fr, "const", targets[["terminal"]]),
                   c_cumulative = {
                       g <- fr[order(fr$yield_cum), ]
                       exp(approx(log(g$yield_cum), log(g$const),
                                  log(targets[["cumulative"]]))$y)
                   })
    }))
})

## --- Phase 2: is the mechanism the one claimed? ----------------------------
message("Phase 2: mechanism controls on eco1")
p1 <- ECOS$eco1$p
alloc_BHP  <- lp_frozen_alloc(p1, "BHP")
alloc_BHPB <- lp_frozen_alloc(p1, "BHPB")

extra <- list(
    BHB          = list(rule = "BHB"),
    frozen_BHP   = list(rule = "frozen", alloc = alloc_BHP),
    frozen_BHPB  = list(rule = "frozen", alloc = alloc_BHPB),
    theta_m0.5   = list(rule = "power", theta = -0.5),
    theta_0.5    = list(rule = "power", theta =  0.5),
    theta_1.5    = list(rule = "power", theta =  1.5))

front2 <- do.call(rbind, mclapply(names(extra), function(nm) {
    a <- extra[[nm]]
    cbind(ecosystem = "eco1",
          cached(paste0("front_eco1_", nm),
                 lp_frontier(p1, a$rule, controls$eco1, NULL,
                             theta = a$theta, alloc = a$alloc,
                             F_ref = 0.1, label = nm)))
}, mc.cores = NCORES))

saveRDS(list(front1 = front1, front2 = front2, alt_calib = alt_calib,
             controls = controls, zfs = zfs),
        "data/robustness/frontiers.rds")
message("wrote data/robustness/frontiers.rds")

## --- Report ----------------------------------------------------------------
all_front <- rbind(front1, front2)

cat("\n== Phase 1: the three published rules, compared at matched terminal yield ==\n")
for (k in names(ECOS)) {
    f <- subset(all_front, ecosystem == k)
    ref <- subset(f, rule == "fixed" & mult == 1)$yield
    cat(sprintf("\n%s  (reference yield %.3f g/m2/yr)\n", k, ref))
    out <- do.call(rbind, lapply(RULES, function(r) {
        fr <- subset(f, rule == r)
        data.frame(rule = r,
                   worst_species  = lp_at_yield(fr, "min_rel_ctrl", ref),
                   n_below_10pct  = lp_at_yield(fr, "n_below_10pct_ctrl", ref),
                   rms_log        = lp_at_yield(fr, "rms_log_ctrl", ref))
    }))
    print(out, row.names = FALSE, digits = 3)
}

cat("\n== Phase 1 item 3: c_P and c_P/B under two calibration criteria ==\n")
print(alt_calib, row.names = FALSE, digits = 3)

cat("\n== Phase 2: mechanism controls on eco1, at the same reference yield ==\n")
ref <- subset(all_front, ecosystem == "eco1" & rule == "fixed" & mult == 1)$yield
out2 <- do.call(rbind, lapply(unique(all_front$rule), function(r) {
    fr <- subset(all_front, ecosystem == "eco1" & rule == r)
    data.frame(rule = r,
               worst_species = lp_at_yield(fr, "min_rel_ctrl", ref),
               n_below_10pct = lp_at_yield(fr, "n_below_10pct_ctrl", ref),
               rms_log       = lp_at_yield(fr, "rms_log_ctrl", ref))
}))
print(out2, row.names = FALSE, digits = 3)
