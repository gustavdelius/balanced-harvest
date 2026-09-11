# docs/robustness.md Phase 5: replication.
#
# Phases 1 and 2 rest on four assembled ecosystems, and between-assemblage
# variation in assembly models is typically large, so a consistent ranking
# across n = 4 is suggestive rather than established.  This assembles twelve
# further ecosystems under exactly the eco1 protocol (no randomisation of the
# search rate, z'_i = 1, F_ref = 0.1) and repeats the frontier comparison on
# all thirteen, reporting distributions rather than exemplars.
#
# Five rules are run: the paper's three, plus the two Phase 2 controls that
# carried the load - F proportional to biomass, and BH_P's year-0 allocation
# frozen for the whole 50 years.
source("R/lp_experiments.R")
source("R/lp_assembly.R")
library(parallel)

NCORES <- as.integer(Sys.getenv("LP_CORES", "4"))
dir.create("data/replication", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/replication", paste0(name, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    val <- force(expr); saveRDS(val, f); val
}

SEEDS <- 301:312
MULTS <- c(0.25, 0.5, 1, 2, 4)

## --- Assemble ---------------------------------------------------------------
message("assembling ", length(SEEDS), " further ecosystems")
new_ecos <- mclapply(SEEDS, function(sd) {
    cached(sprintf("eco_%d", sd), {
        set.seed(sd)
        res <- lp_assemble(dt = 0.01, verbose = FALSE)
        list(params = lp_relax(res$params, years = 50, dt = LP_NUMERICS$dt),
             attempts = res$attempts,
             n_species = nrow(species_params(res$params)))
    })
}, mc.cores = NCORES, mc.preschedule = FALSE)
names(new_ecos) <- paste0("eco", SEEDS)

ecos <- c(list(eco1 = list(params = readRDS("data/ecosystems.rds")$eco1$params,
                           attempts = 23, n_species = 15)),
          new_ecos)

richness <- data.frame(
    ecosystem = names(ecos),
    n_species = sapply(ecos, `[[`, "n_species"),
    attempts  = sapply(ecos, `[[`, "attempts"),
    w_max_min = sapply(ecos, function(e) min(species_params(e$params)$w_max)),
    w_max_max = sapply(ecos, function(e) max(species_params(e$params)$w_max)),
    fished    = sapply(ecos, function(e)
        sum(species_params(e$params)$w_max > LP_FISHING$w_f)))
print(richness, row.names = FALSE, digits = 4)

## --- Frontiers --------------------------------------------------------------
RULESET <- c("fixed", "BHP", "BHPB", "BHB", "frozen_BHP")
spec <- function(nm, params) switch(nm,
    fixed      = list(rule = "fixed"),
    BHP        = list(rule = "BHP"),
    BHPB       = list(rule = "BHPB"),
    BHB        = list(rule = "BHB"),
    frozen_BHP = list(rule = "frozen",
                      alloc = lp_frozen_alloc(params, "BHP")))

message("frontiers: ", length(ecos), " ecosystems x ", length(RULESET), " rules")
jobs <- expand.grid(eco = names(ecos), rule = RULESET, stringsAsFactors = FALSE)
fronts <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
    k <- jobs$eco[i]; nm <- jobs$rule[i]
    p <- ecos[[k]]$params
    ctrl <- cached(paste0("ctrl_", k), lp_control(p))
    a <- spec(nm, p)
    cbind(ecosystem = k,
          cached(sprintf("front_%s_%s", k, nm),
                 lp_frontier(p, a$rule, ctrl, alloc = a$alloc, mults = MULTS,
                             F_ref = 0.1, label = nm)))
}, mc.cores = NCORES, mc.preschedule = FALSE))

## --- Read every rule at its ecosystem's own reference yield -----------------
per_eco <- do.call(rbind, lapply(names(ecos), function(k) {
    f <- subset(fronts, ecosystem == k)
    ref <- subset(f, rule == "fixed" & mult == 1)$yield
    do.call(rbind, lapply(RULESET, function(r) {
        fr <- subset(f, rule == r)
        data.frame(ecosystem = k, rule = r, ref_yield = ref,
                   worst = lp_at_yield(fr, "min_rel_ctrl", ref),
                   rms   = lp_at_yield(fr, "rms_log_ctrl", ref),
                   n10   = lp_at_yield(fr, "n_below_10pct_ctrl", ref))
    }))
}))

saveRDS(list(richness = richness, fronts = fronts, per_eco = per_eco,
             seeds = SEEDS),
        "data/replication/results.rds")

## --- Distributions ----------------------------------------------------------
q <- function(x) sprintf("%.3g [%.3g, %.3g]", median(x, na.rm = TRUE),
                         min(x, na.rm = TRUE), max(x, na.rm = TRUE))
cat("\n== worst-affected species, median [range] over",
    length(ecos), "ecosystems ==\n")
for (r in RULESET)
    cat(sprintf("  %-11s %s\n", r, q(subset(per_eco, rule == r)$worst)))

w <- reshape(per_eco[, c("ecosystem", "rule", "worst")], idvar = "ecosystem",
             timevar = "rule", direction = "wide")
names(w) <- sub("^worst\\.", "", names(w))
cat("\n== BH_P advantage over fixed F (ratio of worst-species outcome) ==\n")
cat(sprintf("  median %.2fx, range %.2f-%.2fx\n",
            median(w$BHP / w$fixed), min(w$BHP / w$fixed), max(w$BHP / w$fixed)))
cat(sprintf("  BH_P beats fixed F in %d of %d ecosystems\n",
            sum(w$BHP > w$fixed), nrow(w)))
cat(sprintf("  BH_P beats BH_P/B in %d of %d\n", sum(w$BHP > w$BHPB), nrow(w)))
cat(sprintf("  full ranking BH_P > fixed > BH_P/B holds in %d of %d\n",
            sum(w$BHP > w$fixed & w$fixed > w$BHPB), nrow(w)))

cat("\n== the two Phase 2 controls, relative to BH_P ==\n")
cat(sprintf("  F ~ B      / BH_P : median %.3f, range %.3f-%.3f\n",
            median(w$BHB / w$BHP), min(w$BHB / w$BHP), max(w$BHB / w$BHP)))
cat(sprintf("  frozen BHP / BH_P : median %.3f, range %.3f-%.3f\n",
            median(w$frozen_BHP / w$BHP), min(w$frozen_BHP / w$BHP),
            max(w$frozen_BHP / w$BHP)))
cat(sprintf("  frozen recovers %.0f%% of the gap on a log scale (median)\n",
            100 * median(log(w$frozen_BHP / w$fixed) / log(w$BHP / w$fixed))))
print(w, row.names = FALSE, digits = 3)
message("wrote data/replication/results.rds")
