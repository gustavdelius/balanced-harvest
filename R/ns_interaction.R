# Build the 16x16 interaction matrix for NS_params plus the rare species.
#
# NS_params' matrix is a spatial-overlap matrix and it is structured, not
# uniform (Sole<->Saithe 0.01, Cod<->Cod 0.79). Rather than give the new
# species a flat constant, each one inherits the overlap profile of the
# resident it most resembles ecologically, optionally scaled by `s`.
suppressMessages(library(mizer))

ANALOGUE <- c("Thornback ray"    = "Gurnard",  # benthic feeder, medium demersal
              "Spurdog"          = "Whiting",  # mobile mid-water piscivore
              "Common skate"     = "Cod",      # large demersal predator
              "Atlantic halibut" = "Saithe")   # large, deeper-water predator

rareInteraction <- function(s = 1, analogue = ANALOGUE) {
    I   <- as.matrix(getInteraction(NS_params))
    res <- rownames(I); new <- names(analogue); all_sp <- c(res, new)
    M <- matrix(NA_real_, length(all_sp), length(all_sp),
                dimnames = list(all_sp, all_sp))
    M[res, res] <- I
    for (x in new) {
        a <- analogue[[x]]
        M[x, res] <- I[a, res]      # x as predator on the residents
        M[res, x] <- I[res, a]      # residents preying on x
    }
    for (x in new) for (y in new)   # among the new species, analogue to analogue
        M[x, y] <- I[analogue[[x]], analogue[[y]]]
    M[, new] <- M[, new] * s        # scale only the new species' exposure/feeding
    M[new, ] <- M[new, ] * s
    M
}
