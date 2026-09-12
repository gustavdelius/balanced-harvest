# Can a fast-slow life-history continuum be assembled?
#
# In the paper's parameterisation the activity factor z_i scales a species'
# search volume and its intrinsic mortality, but not the predation others
# impose on it.  That leaves a steep fitness gradient - d log R0 / d log z
# measured at +3.3 on eco1 - so assembly drives every species to the fastest
# activity on offer and no fast-slow variation survives.  It is the same trap
# that closed the egg-mass axis (see run_egg_assembly.R).
#
# Setting vuln_exp = 1 scales each species' vulnerability to predation with z
# as well, so that EVERY rate it experiences scales together.  Its life is then
# a pure time-rescaling: survivorship is unchanged, reproductive output per
# unit time scales as z, time spent per size interval as 1/z, and R0 is exactly
# invariant (verified: 0.40824 at every z from 0.35 to 2.86).  The axis becomes
# neutral, assembly has nothing to climb, and because P/B = <mu>_B and every
# mortality scales with z, P/B should spread in proportion to the drawn range.
#
# Matched pairs on the same seed:
#   neutral  : vuln_exp = 1, zero fitness gradient
#   gradient : vuln_exp = 0, the paper's treatment, gradient +3.3
#
#   Rscript run_activity_assembly.R <seed> <neutral|gradient>
#
source("R/lp_assembly.R")

args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
mode <- match.arg(args[2], c("neutral", "gradient"))

assembly <- LP_ASSEMBLY
assembly$activity_range <- c(0.35, 2.86)      # ~8-fold, verified neutral
fish <- modifyList(LP_FISH, list(vuln_exp = if (mode == "neutral") 1 else 0))

# The fast end runs ~2.9x the paper's rates, so the time step is cut to match.
DT_ASSEMBLY <- 0.004
DT_RELAX    <- 0.001

dir.create("data/activity", showWarnings = FALSE, recursive = TRUE)
set.seed(seed)
res <- lp_assemble(assembly = assembly, randomise_A = FALSE,
                   dt = DT_ASSEMBLY, fish = fish)
res$params <- lp_relax(res$params, years = 50, dt = DT_RELAX)
res$mode <- mode; res$seed <- seed; res$fish <- fish
saveRDS(res, sprintf("data/activity/%s_%d.rds", mode, seed))

sp <- species_params(res$params)
message(sprintf("%s seed %d: %d species in %d attempts; z %.3f-%.3f",
                mode, seed, nrow(sp), res$attempts, min(sp$z), max(sp$z)))
