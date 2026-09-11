# docs/robustness.md Phase 3 item 8: compensatory recruitment.
#
# The paper deliberately imposes no stock-recruitment relationship: RDD = RDI,
# so recruitment is strictly proportional to egg production and a species whose
# adult mortality exceeds what its lifetime output can sustain declines
# geometrically without limit.  Section 3.2 of the plan argues that this is what
# makes fixed F look catastrophic, and item 8b found the ranking of the rules
# inverting once a brake was restored by contesting the plankton - but there the
# community changed at the same time, so the two could not be separated.
#
# setBevertonHolt() imposes
#
#     RDD = RDI * R_max / (RDI + R_max)
#
# while rescaling erepro so that the *initial state is exactly preserved*.  The
# community, its size structure and its steady state are therefore identical
# across the sweep; only the sensitivity of recruitment to egg production
# changes.  With reproduction_level L = RDD/R_max at the initial state,
#
#     d ln RDD / d ln RDI = R_max / (RDI + R_max) = 1 - L
#
# so L = 0 is the paper (recruitment tracks egg production one for one) and
# L -> 1 is recruitment pinned at R_max regardless of spawning stock.
source("R/lp_experiments.R")
library(parallel)

NCORES <- as.integer(Sys.getenv("LP_CORES", "3"))
dir.create("data/recruitment", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/recruitment", paste0(name, ".rds"))
    if (file.exists(f)) { message("  [cached] ", name); return(readRDS(f)) }
    val <- force(expr); saveRDS(val, f); val
}

LEVELS <- c(0, 0.25, 0.5, 0.75, 0.9)

p0 <- readRDS("data/ecosystems.rds")$eco1$params
zf <- rep(1, nrow(species_params(p0)))

with_level <- function(L) {
    if (L == 0) return(p0)          # noRDD, as published
    setBevertonHolt(p0, reproduction_level = rep(L, nrow(species_params(p0))),
                    info_level = 0)
}

## --- Check that imposing the brake really did leave the state alone --------
message("checking the steady state is preserved")
drift <- do.call(rbind, mclapply(LEVELS, function(L) {
    cached(sprintf("drift_%g", L), {
        p <- with_level(L)
        s <- project(p, t_max = 50, dt = 0.01, t_save = 50,
                     progress_bar = FALSE)
        b <- getBiomass(s)
        data.frame(level = L,
                   erepro = mean(species_params(p)$erepro),
                   recruit_elasticity = 1 - L,
                   max_drift_50yr = max(abs(b[2, ] / b[1, ] - 1)))
    })
}, mc.cores = NCORES))
print(drift, row.names = FALSE, digits = 3)

## --- The three rules under each level of compensation ----------------------
message("harvest comparison at each reproduction level")
jobs <- expand.grid(level = LEVELS, rule = RULES, stringsAsFactors = FALSE)
fronts <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
    L <- jobs$level[i]; r <- jobs$rule[i]
    p <- with_level(L)
    ctrl <- cached(sprintf("ctrl_%g", L), lp_control(p, zf))
    cbind(level = L,
          cached(sprintf("front_%g_%s", L, r),
                 lp_frontier(p, r, ctrl, zf, F_ref = 0.1)))
}, mc.cores = NCORES))

saveRDS(list(drift = drift, fronts = fronts, levels = LEVELS),
        "data/recruitment/results.rds")

cat("\n== three rules at each reproduction level, matched terminal yield ==\n")
summ <- do.call(rbind, lapply(LEVELS, function(L) {
    f <- subset(fronts, level == L)
    ref <- subset(f, rule == "fixed" & mult == 1)$yield
    do.call(rbind, lapply(RULES, function(r) {
        fr <- subset(f, rule == r)
        data.frame(level = L, rule = r, ref_yield = ref,
                   worst_species = lp_at_yield(fr, "min_rel_ctrl", ref),
                   n_below_10pct = lp_at_yield(fr, "n_below_10pct_ctrl", ref),
                   rms_log       = lp_at_yield(fr, "rms_log_ctrl", ref))
    }))
}))
print(summ, row.names = FALSE, digits = 3)
message("wrote data/recruitment/results.rds")
