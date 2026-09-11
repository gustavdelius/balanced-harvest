# Assemble the model ecosystems (Appendix B) and cache them.
#
#   eco1        : the ecosystem used for Figs 2, 3, 4 (no variation in A_i)
#   eco_r1..r3  : three further ecosystems with A_i randomised, for Figs 5, 6
#
# Assembly uses dt = 0.01 yr; the finished ecosystems are then relaxed for a
# further 50 years at the paper's dt = 0.002 yr before any fishing is applied.
source("R/lp_assembly.R")

OUT <- "data/ecosystems.rds"
dir.create("data", showWarnings = FALSE)

assemble_one <- function(seed, randomise_A) {
    set.seed(seed)
    res <- lp_assemble(randomise_A = randomise_A, dt = 0.01)
    message(sprintf("  seed %d: %d species in %d attempts",
                    seed, nrow(species_params(res$params)), res$attempts))
    res$params <- lp_relax(res$params, years = 50, dt = LP_NUMERICS$dt)
    res
}

message("assembling eco1 (Figs 2-4)")
eco1 <- assemble_one(101, randomise_A = FALSE)

ecor <- list()
for (k in 1:3) {
    message("assembling eco_r", k, " (Figs 5, 6)")
    ecor[[k]] <- assemble_one(200 + k, randomise_A = TRUE)
}

saveRDS(list(eco1 = eco1, ecor = ecor), OUT)
message("wrote ", OUT)
