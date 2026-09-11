# Survey biomass for the rare species, from the ICES North Sea IBTS.
#
# Swept-area biomass: catch per haul is divided by the area the net actually
# fished (tow distance x wing spread) and raised to the area of the North Sea.
# Validated against cod, whose ICES stock size is known independently.
#
# An earlier version scaled the rare species relative to NS_params species
# instead. That failed: NS_params gives Dab 10.6 kt and Gurnard 62 kt, while
# the survey catches 660x more dab than gurnard, so the estimates spanned four
# orders of magnitude depending on which reference was chosen.
suppressMessages({library(icesDatras); library(dplyr); library(mizer)})

YEARS <- 2015:2019
APHIA <- c("Thornback ray" = 105883, "Common skate" = 105869,
           "Atlantic halibut" = 127138, "Spurdog" = 105923,
           # references, all present in NS_params
           "Dab" = 127139, "Whiting" = 126438, "Cod" = 126436, "Gurnard" = 127203)
REFS  <- c("Dab", "Whiting", "Cod", "Gurnard")

CACHE <- "data/datras_ns_ibts_q1.rds"
if (file.exists(CACHE)) {
    raw <- readRDS(CACHE)
} else {
    raw <- list(
        hl = bind_rows(lapply(YEARS, function(y)
                 getHLdata("NS-IBTS", y, 1) |> mutate(Year = y))),
        hh = bind_rows(lapply(YEARS, function(y)
                 getHHdata("NS-IBTS", y, 1) |> mutate(Year = y))))
    saveRDS(raw, CACHE)
}

haul_key <- function(d)
    paste(d$Year, d$Country, d$Ship, d$Gear, d$StNo, d$HaulNo, sep = "|")

# Valid hauls only, with their durations
hh <- raw$hh |> filter(HaulVal == "V", HaulDur > 0) |>
      mutate(key = haul_key(cur_data_all())) |> distinct(key, .keep_all = TRUE)

# Length-weight parameters from FishBase for every species used
SCI <- c("Thornback ray" = "Raja clavata", "Common skate" = "Dipturus batis",
         "Atlantic halibut" = "Hippoglossus hippoglossus",
         "Spurdog" = "Squalus acanthias", "Dab" = "Limanda limanda",
         "Whiting" = "Merlangius merlangus", "Cod" = "Gadus morhua",
         "Gurnard" = "Eutrigla gurnardus")
lw <- rfishbase::estimate(unname(SCI), fields = c("Species", "a", "b"))
LW <- lw[match(SCI, lw$Species), c("a", "b")]; rownames(LW) <- names(SCI)

## Biomass caught per hour of trawling, by species -------------------------
cpue <- raw$hl |>
    filter(Valid_Aphia %in% APHIA, SpecVal %in% c(1, 4, 7, 10),
           !is.na(HLNoAtLngt), HLNoAtLngt > 0) |>
    mutate(key     = haul_key(cur_data_all()),
           species = names(APHIA)[match(Valid_Aphia, APHIA)],
           len_cm  = ifelse(LngtCode == "1", LngtClass, LngtClass / 10),
           number  = HLNoAtLngt * ifelse(is.na(SubFactor), 1, SubFactor),
           wt_g    = LW$a[match(species, rownames(LW))] *
                     len_cm ^ LW$b[match(species, rownames(LW))]) |>
    filter(key %in% hh$key) |>
    group_by(species) |>
    summarise(kg = sum(number * wt_g, na.rm = TRUE) / 1000,
              n  = sum(number, na.rm = TRUE), .groups = "drop") |>
    mutate(kg_per_hr = kg / sum(hh$HaulDur / 60))


## Swept area -------------------------------------------------------------
NS_AREA_KM2 <- 570000   # ICES Subarea 4

hh <- hh |>
    mutate(dist = suppressWarnings(as.numeric(Distance)),
           wing = suppressWarnings(as.numeric(WingSpread)),
           dist = ifelse(dist > 0, dist, NA), wing = ifelse(wing > 0, wing, NA),
           # GOV wing spread is tightly set by gear design; impute the median
           wing = ifelse(is.na(wing), median(wing, na.rm = TRUE), wing),
           # where tow distance is missing, reconstruct from speed and duration
           dist = ifelse(is.na(dist),
                         suppressWarnings(as.numeric(GroundSpeed)) * 1852 / 60 * HaulDur,
                         dist),
           swept_km2 = dist * wing / 1e6) |>
    filter(!is.na(swept_km2), swept_km2 > 0)

total_swept <- sum(hh$swept_km2)

cpue <- raw$hl |>
    filter(Valid_Aphia %in% APHIA, SpecVal %in% c(1, 4, 7, 10),
           !is.na(HLNoAtLngt), HLNoAtLngt > 0) |>
    mutate(key     = haul_key(pick(everything())),
           species = names(APHIA)[match(Valid_Aphia, APHIA)],
           len_cm  = ifelse(LngtCode == "1", LngtClass, LngtClass / 10),
           number  = HLNoAtLngt * ifelse(is.na(SubFactor), 1, SubFactor),
           wt_g    = LW$a[match(species, rownames(LW))] *
                     len_cm ^ LW$b[match(species, rownames(LW))]) |>
    filter(key %in% hh$key) |>
    group_by(species) |>
    summarise(kg = sum(number * wt_g, na.rm = TRUE) / 1000,
              n  = sum(number, na.rm = TRUE), .groups = "drop") |>
    mutate(kg_per_km2 = kg / total_swept,
           tonnes     = kg_per_km2 * NS_AREA_KM2 / 1000)

cat("\nNS-IBTS Q1", min(YEARS), "-", max(YEARS), ":", nrow(hh), "hauls,",
    signif(total_swept, 3), "km2 swept (", signif(100*total_swept/NS_AREA_KM2, 2),
    "% of the North Sea )\n\n")
print(cpue |> arrange(desc(tonnes)) |>
      mutate(across(where(is.numeric), ~signif(., 3))) |> as.data.frame(),
      row.names = FALSE)

cat("\nValidation - swept-area estimate vs NS_params model biomass (tonnes):\n")
mb <- getBiomass(NS_params) / 1e6
chk <- cpue |> filter(species %in% names(mb)) |>
    mutate(model = signif(mb[species], 3), swept = signif(tonnes, 3),
           ratio = signif(tonnes / mb[species], 2)) |> select(species, swept, model, ratio)
print(as.data.frame(chk), row.names = FALSE)

obs <- setNames(cpue$tonnes * 1e6, cpue$species)   # grams
saveRDS(list(cpue = cpue, biomass_observed = obs, years = YEARS,
             n_hauls = nrow(hh), swept_km2 = total_swept,
             area_km2 = NS_AREA_KM2), "data/ns_biomass_observed.rds")
