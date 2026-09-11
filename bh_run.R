# Run the three harvesting regimes and cache everything the figures need.
source("bh_model.R")
library(dplyr); library(tidyr)

OUT <- "bh_results.rds"

p0 <- projectToSteady(setup(), t_max = 200, progress_bar = FALSE)  # unfished state
sp <- species_params(p0)$species
set.seed(42)
z <- setNames(runif(length(sp), 0.5, 1.5), sp)   # implementation factor, U(0.5, 1.5)

run <- function(rule, const)
    project(setRule(p0, rule, const, z), t_max = T_MAX, t_save = 1,
            progress_bar = FALSE)
totalYield <- function(sim) sum(getYield(sim)[as.character(T_MAX), ])

## Calibrate c_P and c_P/B to the same total yield as the fixed-F regime -----
sims   <- list(fixed = run("fixed", F_FIXED))
target <- totalYield(sims$fixed)
calibrate <- function(rule, interval)
    exp(uniroot(function(l) totalYield(run(rule, exp(l))) - target,
                interval = log(interval), tol = 1e-3)$root)

const <- c(fixed = F_FIXED,
           BHP   = calibrate("BHP",  c(3e-13, 1e-11)),
           BHPB  = calibrate("BHPB", c(0.05, 3)))
sims$BHP  <- run("BHP",  const[["BHP"]])
sims$BHPB <- run("BHPB", const[["BHPB"]])

message("calibrated: c_P = ", signif(const[["BHP"]], 4),
        ", c_P/B = ", signif(const[["BHPB"]], 4),
        "; target yield = ", signif(target, 4))

## Tidy up the results ------------------------------------------------------
labels <- c(fixed = "Fixed F", BHP = "BH[P]", BHPB = "BH[P/B]")
tidy <- function(m, value)
    as.data.frame.table(m, responseName = value, stringsAsFactors = FALSE) |>
        rename(time = 1, species = 2) |> mutate(time = as.numeric(time))

fished <- sp[production(p0, initialN(p0), getEGrowth(p0)) > 0]

biomass <- bind_rows(lapply(names(sims), function(r)
    tidy(getBiomass(sims[[r]]), "biomass") |> mutate(regime = r)))
fmort <- bind_rows(lapply(names(sims), function(r) {
    m <- FOverTime(sims[[r]]); dimnames(m) <- list(time = 0:T_MAX, species = sp)
    tidy(m, "F") |> mutate(regime = r)
}))
yield <- bind_rows(lapply(names(sims), function(r) {
    s <- sims[[r]]; n <- N(s)[as.character(T_MAX), , ]
    data.frame(species = sp, regime = r,
               yield = getYield(s)[as.character(T_MAX), ],
               B = biomass_fished(p0, n), P = production(p0, n, getEGrowth(p0, n = n)))
}))
spectrum <- function(n, regime)
    as.data.frame.table(sweep(n, 2, w(p0)^2, "*"), responseName = "b",
                        stringsAsFactors = FALSE) |>
        rename(species = 1, w = 2) |>
        mutate(w = as.numeric(w), regime = regime)

spectra <- bind_rows(
    spectrum(initialN(p0), "unfished"),
    lapply(names(sims), function(r)
        spectrum(N(sims[[r]])[as.character(T_MAX), , ], r)))

## Intensity sweep: the yield / depletion trade-off ------------------------
# For each regime, scale the constant over a range and record total yield
# against the most-depleted species, as a measure of risk to biodiversity.
unfished <- getBiomass(sims$fixed)[1, ]
sweep_one <- function(rule, const) {
    sim <- run(rule, const)
    b   <- getBiomass(sim)[as.character(T_MAX), ] / unfished
    data.frame(regime = rule, const = const, yield = totalYield(sim),
               species = names(b), rel = as.numeric(b), row.names = NULL)
}
mult   <- c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6)
sweeps <- bind_rows(lapply(names(const), function(r)
    bind_rows(lapply(mult, function(m) sweep_one(r, const[[r]] * m)))))
message("sweep done")

saveRDS(list(biomass = biomass, fmort = fmort, yield = yield, spectra = spectra, sweeps = sweeps,
             const = const, z = z, target = target, fished = fished, P0 = production(p0, initialN(p0), getEGrowth(p0)),
             B0 = biomass_fished(p0, initialN(p0)),
             species = sp, w_f = W_F, w_max = setNames(species_params(p0)$w_max, sp), t_max = T_MAX,
             unfished_biomass = getBiomass(sims$fixed)[1, ], labels = labels),
        OUT)
message("wrote ", OUT)
