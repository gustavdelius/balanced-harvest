# Does egg mass, made an independent life-history axis, spread P/B across
# species?  See docs/index.md: in the paper's parameterisation P/B is nearly
# species-independent, which is what makes production proportional to biomass
# and what makes BH_P/B almost indistinguishable from a fixed F.
#
# Egg mass is the strongest available lever because sub-gram fish hold ~1% of
# the biomass but generate 18-59% of the production, and the larval mortality
# of Eq. (A.6) acts on absolute size: a species with a 0.1 g egg starts at the
# top of the gauntlet, one with a 1 mg egg runs the whole of it.
#
# Matched pairs: the same seed is assembled with eggs varied and with the
# paper's common w_0, so that any change in the spread of P/B is not just a
# different assemblage.  Usage:
#
#   Rscript run_egg_assembly.R <seed> <varied|fixed>
#
source("R/lp_assembly.R")

args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
mode <- match.arg(args[2], c("varied", "fixed"))

assembly <- LP_ASSEMBLY
# Upwards from the paper's w_0 = 1e-3 g to w_L = 0.1 g, the size at which
# larval mortality dies away.  Upwards because mizer takes the size grid's
# floor from the smallest egg, so going below w_0 would move the grid.
if (mode == "varied") assembly$w_egg_range <- c(LP_FISH$w_egg, LP_FISH$w_L)

dir.create("data/egg", showWarnings = FALSE, recursive = TRUE)
out <- sprintf("data/egg/%s_%d.rds", mode, seed)

set.seed(seed)
res <- lp_assemble(assembly = assembly, randomise_A = FALSE, dt = 0.01)
res$params <- lp_relax(res$params, years = 50, dt = LP_NUMERICS$dt)
res$mode <- mode; res$seed <- seed
saveRDS(res, out)

sp <- species_params(res$params)
message(sprintf("%s seed %d: %d species in %d attempts; eggs %.4g-%.4g g",
                mode, seed, nrow(sp), res$attempts,
                min(sp$w_min), max(sp$w_min)))
