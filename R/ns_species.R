# Pull life-history parameters for the rare species from FishBase via rfishbase
# and assemble a mizer species_params table. Every value is recorded with its
# FishBase source field so it can be checked.
#
# Species: three elasmobranchs that a North Sea mixed demersal fishery catches
# as bycatch and that are genuinely depleted, plus one large slow teleost.
suppressMessages(library(rfishbase)); suppressMessages(library(dplyr))

SPECIES <- c("Raja clavata"              = "Thornback ray",
             "Squalus acanthias"         = "Spurdog",
             "Dipturus batis"            = "Common skate",
             "Hippoglossus hippoglossus" = "Atlantic halibut")
sci <- names(SPECIES)

lw  <- estimate(sci, fields = c("Species", "a", "b"))
mat <- maturity(sci)  |> filter(!is.na(Lm)) |> group_by(Species) |>
       summarise(l_mat = median(Lm, na.rm = TRUE), n_mat = n(), .groups = "drop")
pop <- popgrowth(sci) |> filter(!is.na(Loo)) |> group_by(Species) |>
       summarise(l_max = median(Loo, na.rm = TRUE),
                 k_vb  = median(K, na.rm = TRUE), n_pop = n(), .groups = "drop")

RARE <- data.frame(Species = sci, species = unname(SPECIES)) |>
    left_join(lw, by = "Species") |>
    left_join(mat, by = "Species") |>
    left_join(pop, by = "Species")

print(RARE)
saveRDS(RARE, "data/ns_species_params_raw.rds")

## Assemble the mizer species_params table ----------------------------------
# l_max here is FishBase's asymptotic length Loo, which is what mizer's w_max
# wants. No w_min: addSpecies() cannot take one (sizespectrum/mizer#610), so
# these species are given the model's 1 mg egg rather than live-borne pups.
RARE <- RARE |>
    transmute(species, Scientific_name = Species, a, b, l_max, l_mat, k_vb,
              w_max = a * l_max^b, w_mat = a * l_mat^b)
print(RARE)
saveRDS(RARE, "data/ns_species_params.rds")
