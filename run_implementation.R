# docs/robustness.md Phase 4: the rule as it would actually be implemented.
#
# Phase 5 established that the adaptive feedback, not the initial allocation,
# is what makes BH_P work.  That makes it important how well the feedback
# survives being applied the way a real fishery would apply it: updated every
# few years from survey estimates rather than continuously from perfect
# knowledge.  Two degradations are imposed together,
#
#   interval  years between recalculations of F_i (0 = continuous, the ideal;
#             50 = set once at the start and never revised, the frozen control)
#   sigma     log-scale s.d. of mean-preserving log-normal error on the
#             observed allocation, redrawn independently at every update
#
# on five ecosystems, with a short frontier at each setting so that everything
# is still read at the ecosystem's own reference yield.
source("R/lp_experiments.R")
library(parallel)

NCORES <- as.integer(Sys.getenv("LP_CORES", "4"))
dir.create("data/implementation", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/implementation", paste0(name, ".rds"))
    if (file.exists(f)) return(readRDS(f))
    val <- force(expr); saveRDS(val, f); val
}

ECOS <- c("eco1", paste0("eco", 301:304))
MULTS <- c(0.5, 1, 2, 4)
CONFIGS <- data.frame(
    interval = c(0,  5,  10, 25, 50,  5,   10,  5),
    sigma    = c(0,  0,  0,  0,  0,   0.5, 0.5, 1.0))
CONFIGS$label <- with(CONFIGS, sprintf("k%g_s%g", interval, sigma))

get_params <- function(k) {
    if (k == "eco1") readRDS("data/ecosystems.rds")$eco1$params
    else readRDS(sprintf("data/replication/eco_%s.rds", sub("^eco", "", k)))$params
}
rep_fronts <- readRDS("data/replication/results.rds")$fronts

jobs <- expand.grid(eco = ECOS, cfg = seq_len(nrow(CONFIGS)),
                    stringsAsFactors = FALSE)
message("Phase 4: ", nrow(jobs), " configurations x ", length(MULTS),
        " intensities")

res <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
    k <- jobs$eco[i]; ci <- jobs$cfg[i]
    cfg <- CONFIGS[ci, ]
    p <- get_params(k)
    ctrl <- readRDS(sprintf("data/replication/ctrl_%s.rds", k))
    base <- lp_base_const(p, "BHP", F_ref = 0.1)
    cached(sprintf("%s_%s", k, cfg$label), {
        do.call(rbind, lapply(MULTS, function(m) {
            tr <- lp_harvest_periodic(p, "BHP", base * m, interval = cfg$interval,
                                      sigma = cfg$sigma,
                                      seed = 1000 * ci + which(ECOS == k))
            cbind(data.frame(ecosystem = k, interval = cfg$interval,
                             sigma = cfg$sigma, rule = "BHP", mult = m),
                  lp_metrics(tr, ctrl))
        }))
    })
}, mc.cores = NCORES, mc.preschedule = FALSE))

## --- read every configuration at its ecosystem's reference yield ------------
summ <- do.call(rbind, lapply(ECOS, function(k) {
    f <- subset(rep_fronts, ecosystem == k)
    ref <- subset(f, rule == "fixed" & mult == 1)$yield
    fixed_worst <- lp_at_yield(subset(f, rule == "fixed"), "min_rel_ctrl", ref)
    do.call(rbind, lapply(seq_len(nrow(CONFIGS)), function(ci) {
        cfg <- CONFIGS[ci, ]
        fr <- subset(res, ecosystem == k & interval == cfg$interval &
                          sigma == cfg$sigma)
        w <- lp_at_yield(fr, "min_rel_ctrl", ref)
        data.frame(ecosystem = k, interval = cfg$interval, sigma = cfg$sigma,
                   worst = w, fixed_worst = fixed_worst,
                   advantage = w / fixed_worst)
    }))
}))
saveRDS(list(res = res, summ = summ, configs = CONFIGS),
        "data/implementation/results.rds")

cat("\n== BH_P advantage over fixed F, median [range] over 5 ecosystems ==\n")
cat(sprintf("%9s %7s   %s\n", "interval", "sigma", "advantage"))
for (ci in seq_len(nrow(CONFIGS))) {
    cfg <- CONFIGS[ci, ]
    v <- subset(summ, interval == cfg$interval & sigma == cfg$sigma)$advantage
    cat(sprintf("%7g yr %7g   %.2fx  [%.2f, %.2f]   (n=%d)\n",
                cfg$interval, cfg$sigma, median(v, na.rm = TRUE),
                min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(!is.na(v))))
}
cat("\n== worst-affected species by configuration ==\n")
print(reshape(summ[, c("ecosystem","interval","sigma","worst")],
              idvar = c("interval","sigma"), timevar = "ecosystem",
              direction = "wide"), row.names = FALSE, digits = 3)
message("wrote data/implementation/results.rds")
