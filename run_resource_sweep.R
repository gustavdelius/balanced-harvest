# docs/robustness.md section 3.6 / Phase 3 item 8b.
#
# Does switching on the larval competition that Appendix B nominates as the
# model's density-dependent brake change the comparison between harvest rules?
#
# Primary production is held fixed (r_0 * a_0 = 20000) so that the variants
# differ in how *contested* the plankton is, not in how much of it there is.
# Ecosystems are assembled from the same seed as eco1, so the invader draws are
# identical across variants and only the plankton dynamics differ.
source("R/lp_resource_variants.R")
source("R/lp_experiments.R")
library(parallel)

NCORES <- as.integer(Sys.getenv("LP_CORES", "2"))
dir.create("data/resource", showWarnings = FALSE, recursive = TRUE)
cached <- function(name, expr) {
    f <- file.path("data/resource", paste0(name, ".rds"))
    if (file.exists(f)) { message("  [cached] ", name); return(readRDS(f)) }
    val <- force(expr); saveRDS(val, f); val
}

# r_0 = 10 is the published value; 3 is the fastest slowdown that still has a
# stable fixed point in the single-species screen; 1 is past the Hopf.
VARIANTS <- list(published = 10, slow3 = 3, slow10 = 1)

## --- Assemble one ecosystem per variant ------------------------------------
message("assembling ecosystems (same seed as eco1, so the draws match)")
ecos <- mclapply(names(VARIANTS), function(nm) {
    cached(paste0("eco_", nm), {
        r0 <- VARIANTS[[nm]]
        pl <- lp_plankton_variant(r0, 2e4 / r0, 2e3)
        set.seed(101)
        res <- lp_assemble(plankton = pl, dt = 0.01, verbose = FALSE)
        list(params = lp_relax(res$params, years = 50, dt = LP_NUMERICS$dt),
             attempts = res$attempts, n_species = nrow(species_params(res$params)))
    })
}, mc.cores = min(NCORES, 3))
names(ecos) <- names(VARIANTS)

## --- Does the assembled system settle, and is larval competition on? -------
message("diagnosing")
diag <- do.call(rbind, mclapply(names(VARIANTS), function(nm) {
    cached(paste0("diag_", nm), {
        p <- ecos[[nm]]$params
        # 200 unfished years: does total biomass settle or cycle?
        s <- project(p, t_max = 200, dt = 0.01, t_save = 5, progress_bar = FALSE)
        tb <- rowSums(getBiomass(s))
        late <- tb[as.numeric(names(tb)) >= 100]
        d <- lp_resource_diagnostics(p)
        data.frame(variant = nm, r0 = VARIANTS[[nm]],
                   n_species = ecos[[nm]]$n_species,
                   attempts = ecos[[nm]]$attempts,
                   as.list(d),
                   cv_late = sd(late) / mean(late),
                   swing_late = max(late) / min(late))
    })
}, mc.cores = min(NCORES, 3)))

print(diag, row.names = FALSE, digits = 3)

## --- The three-rule comparison under each variant --------------------------
# Contesting the plankton excludes large-bodied species entirely (see the
# diagnostics above), so under the slowed variants no species reaches the
# paper's w_f = 400 g and there is no fishery at all.  To compare the rules
# the fishery has to enter at the same point in a typical species' life
# history in each community, so w_f is scaled to hold w_f / median(w_max)
# at the value it takes in the published ecosystem.
# Matching w_f on body mass is meaningless here: with the plankton contested,
# juvenile growth is crushed and the top third of community biomass sits at
# 12 g (r0 = 3) or 1.6 g (r0 = 1) rather than 542 g, so a 400 g fishery would
# catch nothing at all.  Instead w_f is set so that the fishery is exposed to
# the same FRACTION of community biomass in every variant as it is in the
# published one.
WF_FRAC <- lp_fraction_above(ecos$published$params, LP_FISHING$w_f)
w_fs <- sapply(ecos, function(e) lp_wf_for_fraction(e$params, WF_FRAC))
message(sprintf("matching exposed biomass fraction %.3f", WF_FRAC))
message("w_f per variant: ",
        paste(sprintf("%s = %.0f g", names(w_fs), w_fs), collapse = ", "))

message("harvest comparison under each variant")
MULTS <- c(0.25, 0.5, 1, 2, 4)
jobs <- expand.grid(variant = names(VARIANTS), rule = RULES,
                    stringsAsFactors = FALSE)
fronts <- do.call(rbind, mclapply(seq_len(nrow(jobs)), function(i) {
    nm <- jobs$variant[i]; r <- jobs$rule[i]
    p <- ecos[[nm]]$params
    wf <- w_fs[[nm]]
    ctrl <- cached(paste0("ctrl_", nm), lp_control(p, w_f = wf))
    cbind(variant = nm, w_f = wf,
          cached(sprintf("front_%s_%s", nm, r),
                 lp_frontier(p, r, ctrl, mults = MULTS, F_ref = 0.1,
                             w_f = wf)))
}, mc.cores = NCORES))

saveRDS(list(diag = diag, fronts = fronts, variants = VARIANTS, w_f = w_fs),
        "data/resource/results.rds")

cat("\n== three rules under each plankton variant, at matched terminal yield ==\n")
for (nm in names(VARIANTS)) {
    f <- subset(fronts, variant == nm)
    ref <- subset(f, rule == "fixed" & mult == 1)$yield
    cat(sprintf("\n%s (r0 = %g, w_f = %.0f g), reference yield %.4f g/m2/yr\n",
                nm, VARIANTS[[nm]], w_fs[[nm]], ref))
    print(do.call(rbind, lapply(RULES, function(r) {
        fr <- subset(f, rule == r)
        data.frame(rule = r,
                   worst_species = lp_at_yield(fr, "min_rel_ctrl", ref),
                   n_below_10pct = lp_at_yield(fr, "n_below_10pct_ctrl", ref),
                   rms_log       = lp_at_yield(fr, "rms_log_ctrl", ref))
    })), row.names = FALSE, digits = 3)
}
message("wrote data/resource/results.rds")
