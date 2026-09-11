# ---------------------------------------------------------------------------
# Figures for the North Sea exercise. Reads data/ns_rare.rds and
# data/ns_biomass_observed.rds; writes PNGs into docs/figures/.
# ---------------------------------------------------------------------------
library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork)

RULE_LAB <- c(none = "no fishing", fixed = "F[i]*': constant'",
              BHP  = "F[i]*': '*BH[P]", BHPB = "F[i]*': '*BH[P/B]")
RULE_COL <- c(none = "#8a8984", fixed = "#2a78d6", BHP = "#eb6834", BHPB = "#1baf7a")

theme_ns <- function(base = 9) {
    theme_bw(base_size = base) +
        theme(panel.grid = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = base, hjust = 0.5),
              legend.position = "top", legend.title = element_blank(),
              legend.key.height = unit(9, "pt"))
}
asRule <- function(d) mutate(d, regime = factor(regime, levels = names(RULE_COL)))

r   <- readRDS("data/ns_rare.rds")
obs <- readRDS("data/ns_biomass_observed.rds")
dir.create("docs/figures", showWarnings = FALSE, recursive = TRUE)

## Fig 1 - how each rule allocates F across the community ------------------
# Species whose maximum size is below the 400 g entry size are never actually
# caught; a constant-F rule still assigns them a nominal rate, so drop them.
fished <- r$sweep |> filter(regime == "BHP", F50 > 0) |> pull(species) |> unique()
f50 <- r$sweep |> filter(rl == 0.05, regime != "none", species %in% fished) |>
    mutate(grp = ifelse(rare, "rare (survey-derived)", "NS_params resident")) |>
    asRule()
ord <- f50 |> filter(regime == "fixed") |> arrange(F50) |> pull(species)

fig1 <- ggplot(f50, aes(F50, factor(species, levels = ord),
                        colour = regime, shape = grp)) +
    geom_line(aes(group = species), colour = "grey80", linewidth = 1.6,
              lineend = "round") +
    geom_point(size = 1.9) +
    scale_colour_manual(values = RULE_COL, labels = parse(text = RULE_LAB[-1]),
                        breaks = names(RULE_COL)[-1]) +
    scale_shape_manual(values = c(16, 17)) +
    scale_x_log10(labels = label_log()) +
    labs(x = expression(italic(F)[i]~at~year~50~(yr^-1)), y = NULL) +
    theme_ns() + guides(colour = guide_legend(order = 1))
ggsave("docs/figures/ns_fig1.png", fig1, width = 6.4, height = 4.0, dpi = 200)

## Fig 2 - the rare species, against the no-fishing control ----------------
rel <- r$sweep |> filter(rare) |> select(rl, regime, species, rel) |>
    pivot_wider(names_from = regime, values_from = rel) |>
    mutate(across(c(fixed, BHP, BHPB), ~ . / none)) |>
    select(-none) |>
    pivot_longer(c(fixed, BHP, BHPB), names_to = "regime", values_to = "rel") |>
    asRule()

fig2 <- ggplot(rel, aes(rl, rel, colour = regime)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    facet_wrap(~species) +
    scale_colour_manual(values = RULE_COL, labels = parse(text = RULE_LAB[-1]),
                        breaks = names(RULE_COL)[-1]) +
    scale_x_log10(labels = label_log()) +
    scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
    labs(x = "reproduction level of the rare species",
         y = "year-50 biomass, relative to\nthe no-fishing control") +
    theme_ns()
ggsave("docs/figures/ns_fig2.png", fig2, width = 6.4, height = 3.4, dpi = 200)

## Fig 3 - the trade-off ---------------------------------------------------
tr <- r$sweep |> filter(rl == 0.05) |> group_by(regime) |>
    summarise(residents = min(rel[!rare]),
              rare = min(rel[rare] / rel[rare][1]), .groups = "drop")
ctrl <- r$sweep |> filter(rl == 0.05, regime == "none") |>
    summarise(res = min(rel[!rare]), rare = min(rel[rare])) 
tr <- r$sweep |> filter(rl == 0.05) |>
    left_join(r$sweep |> filter(rl == 0.05, regime == "none") |>
                  select(species, ctrl = rel), by = "species") |>
    filter(regime != "none") |>
    group_by(regime) |>
    summarise(residents = min(rel[!rare] / ctrl[!rare]),
              rare      = min(rel[rare]  / ctrl[rare]),
              yield     = sum(yield), .groups = "drop") |> asRule()

fig3 <- ggplot(tr, aes(residents, rare, colour = regime)) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 1, linetype = "22", colour = "grey60", linewidth = 0.3) +
    geom_point(size = 3.2) +
    ggrepel::geom_text_repel(aes(label = c(fixed = "constant F", BHP = "BH_P",
                                           BHPB = "BH_P/B")[as.character(regime)]),
                             size = 2.9, show.legend = FALSE, seed = 1) +
    scale_colour_manual(values = RULE_COL, guide = "none") +
    scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
    labs(x = "worst-affected resident, relative to no-fishing control",
         y = "worst-affected rare species,\nrelative to no-fishing control") +
    theme_ns()
ggsave("docs/figures/ns_fig3.png", fig3, width = 5.2, height = 3.6, dpi = 200)

message("wrote docs/figures/ns_fig1-3.png")
