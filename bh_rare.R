# Do genuinely rare, low-productivity species change the verdict on BH_P?
#
# Baseline is the PRESENT-DAY FISHED state, not an unfished reconstruction:
# NS_params ships fished and calibrated to that state, so it is the one
# reference actually fitted to observations. Each scenario asks what happens
# if the fishery switches to rule X for 50 years from now.
#
# Rare species enter with the analogy-based interaction matrix (bh_interaction.R),
# not addSpecies()' flat default of 1, which would expose them more strongly
# than any fitted pair in the model and makes them spuriously non-viable.
#
# Spurdog is excluded: viviparous with ~60 g pups, it cannot be represented
# with the model's 1 mg egg (erepro 2100) and needs per-species w_min, which is
# blocked by sizespectrum/mizer#610.

source("bh_model.R")
source("bh_interaction.R")
library(dplyr)

RARE    <- readRDS("rare_species_params.rds")
OBS     <- local({
    o <- readRDS("rare_biomass_observed.rds")   # swept-area, NS-IBTS Q1 2015-2019
    # NS_params sits 4.4x above the swept-area estimate for cod, the best
    # validated species, so observations are put on the model's own scale.
    o$biomass_observed / (o$biomass_observed[["Cod"]] / getBiomass(NS_params)[["Cod"]])
})
# Spurdog: viviparous with ~60 g pups, needs per-species w_min (mizer#610).
# Common skate: at its observed 2500 t it needs erepro 29, i.e. the model says
# it cannot sustain that biomass by local reproduction. Consistent with what
# the survey shows - 22 individuals in 1831 hauls, all 69-115 cm, none near
# its 145 cm maturity - so a remnant rather than a self-sustaining population.
EXCLUDE <- c("Spurdog", "Common skate")
RARE_SP <- setdiff(RARE$species, EXCLUDE)

# reproduction_level -> 0 is near-density-independent reproduction (no
# compensatory recruitment); the twelve residents sit at ~0.99. Above 0.8 the
# common skate needs erepro > 1, so the sweep stops there.
RL <- c(0.001, 0.01, 0.05, 0.1, 0.2, 0.4)

# Residents are first brought to the model's zero-fishing steady state, and the
# rare species are then added AT that state. This matters: addSpecies() tunes
# erepro to hold a species where it is placed, so adding them at NS_params'
# fished state and then changing the fishery makes them drift out of balance
# (the thornback ray goes extinct unfished). The zero-fishing state is used
# because it is the only stationary starting point available to all regimes
# equally -- it is not a claim about a real unfished North Sea.
residents <- local({
    p <- NULL
    function() {
        if (is.null(p)) p <<- projectToSteady(setup(), t_max = 300,
                                              progress_bar = FALSE)
        p
    }
})

build <- function(rl) {
    p <- addSpecies(residents(), RARE, interaction = rareInteraction(1),
                    info_level = 0)
    p <- removeSpecies(p, EXCLUDE)
    sp <- species_params(p); sp$biomass_observed <- NA
    sp$biomass_observed[match(RARE_SP, sp$species)] <- OBS[RARE_SP]
    species_params(p) <- sp
    # Shape first, then level, then the erepro that sustains it. Without the
    # steadySingleSpecies() step matchBiomasses() rescales a spectrum that is
    # not a steady-state shape, and the species drifts even with no fishing.
    p <- steadySingleSpecies(p, species = RARE_SP)
    p <- matchBiomasses(p, species = RARE_SP, info_level = 0)
    reproduction_level(p) <- setNames(rep(rl, length(RARE_SP)), RARE_SP)
    setup(p)
}

runRegime <- function(p, rule, const, z)
    project(setRule(p, rule, const, z), t_max = T_MAX, t_save = 1,
            progress_bar = FALSE)
totalYield <- function(sim) sum(getYield(sim)[as.character(T_MAX), ])

calibrateAll <- function(p, z) {
    target <- totalYield(runRegime(p, "fixed", F_FIXED, z))
    solve1 <- function(rule, interval)
        exp(uniroot(function(l) totalYield(runRegime(p, rule, exp(l), z)) - target,
                    interval = log(interval), tol = 1e-3)$root)
    list(target = target,
         const = c(none = 0, fixed = F_FIXED,
                   BHP   = solve1("BHP",  c(1e-13, 3e-11)),
                   BHPB  = solve1("BHPB", c(0.02, 5))))
}

sweepRare <- function() {
    set.seed(42)
    p_ref <- build(0.2)
    sp    <- species_params(p_ref)$species
    z     <- setNames(runif(length(sp), 0.5, 1.5), sp)
    cal   <- calibrateAll(p_ref, z)
    message("calibrated: c_P = ", signif(cal$const[["BHP"]], 4),
            ", c_P/B = ", signif(cal$const[["BHPB"]], 4),
            "; target = ", signif(cal$target, 4))

    out <- bind_rows(lapply(RL, function(rl) {
        p <- build(rl)
        bind_rows(lapply(names(cal$const), function(rule) {
            sim <- runRegime(p, rule, cal$const[[rule]], z)
            b   <- getBiomass(sim)
            data.frame(rl = rl, regime = rule, species = colnames(b),
                       rel = as.numeric(b[as.character(T_MAX), ] / b["0", ]),
                       F50 = FOverTime(sim)[T_MAX + 1, ],
                       yield = as.numeric(getYield(sim)[as.character(T_MAX), ]),
                       rare = colnames(b) %in% RARE_SP, row.names = NULL)
        }))
    }))
    list(sweep = out, const = cal$const, target = cal$target, z = z,
         rare = RARE_SP, rare_params = RARE, w_f = W_F, t_max = T_MAX, rl = RL)
}

if (sys.nframe() == 0) {
    res <- sweepRare()
    saveRDS(res, "bh_rare_results.rds")
    message("wrote bh_rare_results.rds")
}
