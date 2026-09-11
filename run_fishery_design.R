# docs/robustness.md Phase 3 items 10 and 11: the fishery's design.
#
#   item 10  w_f, the shared knife-edge entry mass.  The paper uses 400 g for
#            every species and flags it as a simplification.  Combined with a
#            fixed w_mat = w_max/10 it means a 1 kg species is fished only as
#            an adult while a 40 kg species is fished for a decade before it
#            breeds, so part of "large species are vulnerable" may be built in.
#
#   item 11  the size range over which P_i and B_i are measured.  The paper
#            sets it equal to the harvested range on data-availability grounds,
#            which makes P_i exclude juvenile production - most of a species'
#            somatic production - and makes it dominated by the boundary influx
#            at w_f for species whose w_max is not far above it.  Here the
#            measurement range is instead the whole life cycle, with the
#            harvested range left at 400 g, so only the rule's information
#            changes, not what is caught.
source("R/lp_experiments.R")
library(parallel)

NCORES <- as.integer(Sys.getenv("LP_CORES", "4"))
dir.create("data/fishery", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/fishery", paste0(name, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    val <- force(expr); saveRDS(val, f); val
}

ECOS  <- c("eco1", paste0("eco", 301:304))
MULTS <- c(0.5, 1, 2)
RULES3 <- c("fixed", "BHP", "BHPB")
# (w_f, measurement range).  NA measurement means "same as harvested", as
# published; 0.001 g means the whole life cycle from egg size up.
CONFIGS <- data.frame(w_f = c(100, 200, 400, 800, 400),
                      w_m = c(NA, NA, NA, NA, 0.001))
CONFIGS$label <- with(CONFIGS, sprintf("wf%g_%s", w_f,
                                       ifelse(is.na(w_m), "harv", "full")))

get_params <- function(k) {
    if (k == "eco1") return(readRDS("data/ecosystems.rds")$eco1$params)
    readRDS(sprintf("data/replication/eco_%s.rds", sub("^eco", "", k)))$params
}

jobs <- expand.grid(eco = ECOS, cfg = seq_len(nrow(CONFIGS)), rule = RULES3,
                    stringsAsFactors = FALSE)
# fixed F ignores P and B, so the full-measurement variant of it is redundant
jobs <- subset(jobs, !(rule == "fixed" & !is.na(CONFIGS$w_m[cfg])))
message("items 10 and 11: ", nrow(jobs), " frontiers")

fronts <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
    k <- jobs$eco[i]; ci <- jobs$cfg[i]; r <- jobs$rule[i]
    cfg <- CONFIGS[ci, ]
    p <- get_params(k)
    ctrl <- readRDS(sprintf("data/replication/ctrl_%s.rds", k))
    cbind(ecosystem = k, config = cfg$label, w_f = cfg$w_f,
          measure = ifelse(is.na(cfg$w_m), "harvested", "whole life cycle"),
          cached(sprintf("%s_%s_%s", k, cfg$label, r),
                 lp_frontier(p, r, ctrl, mults = MULTS, F_ref = 0.1,
                             w_f = cfg$w_f, w_measure = cfg$w_m)))
}, mc.cores = NCORES, mc.preschedule = FALSE))

## --- read each configuration at its own reference yield ---------------------
summ <- do.call(rbind, lapply(ECOS, function(k) {
    do.call(rbind, lapply(seq_len(nrow(CONFIGS)), function(ci) {
        cfg <- CONFIGS[ci, ]
        # fixed F's frontier is shared between the two measurement variants
        fixed_cfg <- if (is.na(cfg$w_m)) cfg$label else sprintf("wf%g_harv", cfg$w_f)
        ff <- subset(fronts, ecosystem == k & config == fixed_cfg & rule == "fixed")
        if (!nrow(ff)) return(NULL)
        ref <- subset(ff, mult == 1)$yield
        do.call(rbind, lapply(RULES3, function(r) {
            fr <- if (r == "fixed") ff else
                subset(fronts, ecosystem == k & config == cfg$label & rule == r)
            if (!nrow(fr)) return(NULL)
            data.frame(ecosystem = k, w_f = cfg$w_f,
                       measure = ifelse(is.na(cfg$w_m), "harvested", "whole"),
                       rule = r, ref_yield = ref,
                       worst = lp_at_yield(fr, "min_rel_ctrl", ref))
        }))
    }))
}))
saveRDS(list(fronts = fronts, summ = summ, configs = CONFIGS),
        "data/fishery/results.rds")

report <- function(d, title) {
    cat("\n==", title, "==\n")
    w <- reshape(d[, c("ecosystem","key","rule","worst")],
                 idvar = c("ecosystem","key"), timevar = "rule",
                 direction = "wide")
    names(w) <- sub("^worst\\.", "", names(w))
    for (kk in unique(w$key)) {
        x <- subset(w, key == kk)
        cat(sprintf("  %-22s BH_P %.3f [%.3f,%.3f]   advantage over fixed %.2fx [%.2f,%.2f]   BH_P wins %d/%d\n",
                    kk, median(x$BHP), min(x$BHP), max(x$BHP),
                    median(x$BHP / x$fixed), min(x$BHP / x$fixed),
                    max(x$BHP / x$fixed), sum(x$BHP > x$fixed), nrow(x)))
    }
}
a <- subset(summ, measure == "harvested"); a$key <- paste0("w_f = ", a$w_f, " g")
report(a, "item 10: entry size")
b <- subset(summ, w_f == 400); b$key <- paste0("P,B measured over ", b$measure)
report(b, "item 11: measurement range (w_f = 400 g)")
message("wrote data/fishery/results.rds")
