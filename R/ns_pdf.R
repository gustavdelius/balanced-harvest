# Figures for the North Sea PDF. Reads data/ns_residents.rds.
library(ggplot2); library(dplyr); library(tidyr); library(scales); library(ggrepel)

r <- readRDS("data/ns_residents.rds")

## Design tokens ------------------------------------------------------------
SURFACE <- "#fcfcfb"; INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8984"
PAL <- c(fixed = "#2a78d6", BHP = "#eb6834", BHPB = "#1baf7a")
LAB <- c(fixed = "Fixed F", BHP = expression(BH[P]), BHPB = expression(BH[P/B]))
LABC <- c(fixed = "Fixed F", BHP = "BH_P", BHPB = "BH_P/B")

theme_bh <- function(base = 10)
    theme_minimal(base_size = base) +
    theme(plot.background  = element_rect(fill = SURFACE, colour = NA),
          panel.background = element_rect(fill = SURFACE, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
          axis.text   = element_text(colour = INK2),
          axis.title  = element_text(colour = INK2),
          strip.text  = element_text(colour = INK, face = "bold", size = base - 1),
          plot.title    = element_text(colour = INK, face = "bold", size = base + 4),
          plot.subtitle = element_text(colour = INK2, size = base, lineheight = 1.2),
          plot.caption  = element_text(colour = MUTED, size = base - 2, hjust = 0),
          legend.position = "top", legend.title = element_blank(),
          legend.text = element_text(colour = INK2),
          plot.margin = margin(14, 16, 12, 14))

regime_scale <- list(
    scale_colour_manual(values = PAL, labels = LAB, breaks = names(PAL)),
    scale_fill_manual(values = PAL, labels = LAB, breaks = names(PAL)))

# Species ordered by unfished somatic production over the harvested range
ord   <- names(sort(r$P0[r$fished], decreasing = TRUE))
asF   <- function(d) mutate(d, regime = factor(regime, levels = names(PAL)))
onlyF <- function(d) filter(d, species %in% r$fished) |>
                     mutate(species = factor(species, levels = ord))

## Fig 1 — the shape of the rule: F is independent of size above w_f --------
f50 <- r$fmort |> filter(time == r$t_max) |> onlyF() |> asF()
grid <- f50 |> rowwise() |>
    mutate(w = list(10^seq(log10(r$w_f), log10(r$w_max[[as.character(species)]]),
                           length.out = 120))) |>
    tidyr::unnest(w) |> ungroup()
ends <- grid |> group_by(regime, species) |> slice_max(w, n = 1) |> ungroup()

fig1 <- ggplot(grid, aes(w, F, group = species)) +
    annotate("rect", xmin = 1, xmax = r$w_f, ymin = -Inf, ymax = Inf,
             fill = "#e4e3dc", alpha = 0.85) +
    annotate("text", x = 1.6, y = min(grid$F) * 1.4, label = "no fishing",
             hjust = 0, size = 2.6, colour = MUTED) +
    geom_line(colour = INK, linewidth = 0.8, alpha = 0.85) +
    geom_point(data = ends, colour = INK, size = 1.1) +
    geom_text_repel(data = ends, aes(label = species), size = 2.7, colour = INK2,
                    direction = "y", hjust = 0, nudge_x = 0.25, xlim = c(log10(5e4), Inf),
                    segment.colour = "#cfcec9", min.segment.length = 0,
                    max.overlaps = Inf) +
    facet_wrap(~regime, labeller = as_labeller(LABC)) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale()),
                  breaks = 10^(0:4), limits = c(1, 1.2e6),
                  expand = expansion(mult = c(0, 0))) +
    scale_y_log10(labels = label_log()) +
    labs(title = "Fishing mortality does not depend on body size",
         subtitle = paste0("Every species enters one mixed fishery at w_f = ", r$w_f,
             " g, and all fish above that size are caught at the same rate.\n",
             "The three rules differ only in how that one rate is set for each species. Dots mark maximum body size."),
         x = "body mass (g)", y = expression(italic(F)[i]~(yr^-1)),
         caption = "North Sea model (mizer NS_params), year 50. Five species with a maximum size below w_f are never caught and are not shown.") +
    theme_bh()

## Fig 2 — how the single rate is allocated across species ------------------
prod <- data.frame(species = factor(ord, levels = rev(ord)),
                   P = r$P0[ord], B = r$B0[ord])
p2a <- ggplot(prod, aes(P, species)) +
    geom_point(colour = INK, size = 2.6) +
    scale_x_log10(labels = label_log()) +
    labs(title = "Somatic production, unfished",
         x = expression(italic(P)[i]~(g~yr^-1)), y = NULL) +
    theme_bh() + theme(plot.title = element_text(size = 10))

p2b <- ggplot(asF(mutate(f50, species = factor(species, levels = rev(ord)))),
              aes(F, species, colour = regime)) +
    geom_line(aes(group = species), colour = "#d8d7d1", linewidth = 2.5,
              lineend = "round") +
    geom_point(size = 2.6) +
    regime_scale +
    scale_x_log10(labels = label_log()) +
    labs(title = "Fishing mortality at year 50",
         x = expression(italic(F)[i]~(yr^-1)), y = NULL) +
    theme_bh() + theme(plot.title = element_text(size = 10))

fig2 <- patchwork::wrap_plots(p2a, p2b, widths = c(1, 1.5)) +
    patchwork::plot_annotation(
        title = "Balancing is across species, not across sizes",
        subtitle = paste0("Species ordered by unfished somatic production over the harvested range. ",
            "BH_P ties F to production, so\nthe least productive species (Gurnard) is fished at ",
            "1e-5 per year - four orders of magnitude below its fixed-F rate."),
        caption = "All three rules are calibrated to the same total yield in year 50.",
        theme = theme_bh())

## Fig 3 — the adaptive feedback -------------------------------------------
fig3 <- ggplot(asF(onlyF(r$fmort)), aes(time, F, colour = regime)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~species, scales = "free_y", nrow = 2) +
    regime_scale +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Under BH_P the fishing mortality is a feedback, not a setting",
         subtitle = paste0("F_i(t) = c_P P_i(t) responds to the stock's own production as fishing changes it. ",
             "Fixed F is flat by\nconstruction; the balanced rules move to a new level over the first decade and then hold it."),
         x = "year", y = expression(italic(F)[i]~(yr^-1)),
         caption = "Note the independent vertical scales.") +
    theme_bh()

## Fig 4 — what it does to the stocks --------------------------------------
rel <- r$biomass |> mutate(rel = biomass / r$unfished_biomass[species]) |>
    mutate(species = factor(species, levels = r$species)) |> asF()
fig4 <- ggplot(rel, aes(time, rel, colour = regime)) +
    geom_hline(yintercept = 1, colour = MUTED, linewidth = 0.3, linetype = "22") +
    geom_line(linewidth = 0.8) +
    facet_wrap(~species, nrow = 3) +
    regime_scale +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       breaks = seq(0.5, 1.3, 0.2)) +
    labs(title = "Stock response over 50 years of fishing, at equal total yield",
         subtitle = paste0("Biomass relative to the unfished state. The first five species are below w_f ",
             "in maximum size and are never\ncaught directly - what moves them is predation released or intensified by fishing on the others."),
         x = "year", y = "biomass, % of unfished",
         caption = "North Sea model (mizer NS_params).") +
    theme_bh()

## Fig 5 — yield against biomass (cf. their Fig. 4b) ------------------------
yb <- r$yield |> filter(species %in% r$fished, yield > 0) |> asF()
ref <- data.frame(B = range(yb$B))
fig5 <- ggplot(yb, aes(B, yield, colour = regime)) +
    geom_line(data = transform(ref, yield = 3e-13 * B^2), aes(B, yield),
              inherit.aes = FALSE, colour = MUTED, linetype = "22", linewidth = 0.4) +
    geom_line(data = transform(ref, yield = 0.22 * B), aes(B, yield),
              inherit.aes = FALSE, colour = MUTED, linetype = "42", linewidth = 0.4) +
    geom_point(size = 2.4) +
    geom_text_repel(data = filter(yb, regime == "BHP"), aes(label = species),
                    size = 2.6, colour = INK2, show.legend = FALSE,
                    min.segment.length = 0, segment.colour = "#cfcec9") +
    annotate("text", x = 5e8, y = 3e-13 * (5e8)^2 * 4.5, label = "slope 2",
             colour = MUTED, size = 2.8) +
    annotate("text", x = 5e8, y = 0.22 * 5e8 * 2.2, label = "slope 1",
             colour = MUTED, size = 2.8) +
    regime_scale +
    scale_x_log10(labels = label_log()) + scale_y_log10(labels = label_log()) +
    labs(title = "Why BH_P takes so little from the unproductive species",
         subtitle = paste0("Yield against harvestable biomass in year 50. BH_P makes Y_i = c_P P_i B_i, ",
             "and since production and biomass\nrun together this puts the species near a line of slope 2: ",
             "halving a species' biomass quarters its catch."),
         x = expression(harvestable~biomass~italic(B)[i]~(g)),
         y = expression(yield~italic(Y)[i]~(g~yr^-1)),
         caption = "Reference lines are guides, not fits. A fixed F gives slope 1 by construction.") +
    theme_bh()

## Fig 6 — the community size spectrum -------------------------------------
spec <- r$spectra |> group_by(regime, w) |> summarise(b = sum(b), .groups = "drop") |>
    mutate(regime = factor(regime, levels = c("unfished", names(PAL))))
fig6 <- ggplot(spec, aes(w, b, colour = regime)) +
    geom_vline(xintercept = r$w_f, linetype = "22", colour = MUTED, linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c(unfished = INK, PAL),
                        labels = c(unfished = "Unfished", LAB),
                        breaks = c("unfished", names(PAL))) +
    scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
    scale_y_log10(labels = label_log()) +
    coord_cartesian(xlim = c(1, 4e4), ylim = c(2.5e10, 6e12)) +
    labs(title = "Fishing the big fish lifts the middle of the spectrum",
         subtitle = paste0("Total fish biomass density across all twelve species, year 50. The unfished state piles biomass up at maximum\n",
             "body size, where growth ceases and only background mortality applies. Removing that pile cuts predation on\n",
             "1-10 kg fish by more than the fishing it adds - Cod's death rate at 5.7 kg falls from 0.54 to 0.13/yr while F rises\n",
             "by 0.20 - so the middle of the spectrum rises even though total biomass falls."),
         x = "body mass (g)", y = expression(biomass~density~italic(N)(w)~w^2),
         caption = "Dashed line marks w_f = 400 g. Recruitment is essentially unchanged in all three regimes.") +
    theme_bh()

## Fig 7 — the trade-off, swept over fishing intensity ----------------------
sw <- r$sweeps |> group_by(regime, const, yield) |>
    summarise(`most depleted species` = min(rel),
              `Gurnard (least productive)` = rel[species == "Gurnard"],
              .groups = "drop") |>
    pivot_longer(-c(regime, const, yield), names_to = "metric", values_to = "rel") |>
    mutate(metric = factor(metric, levels = c("most depleted species",
                                              "Gurnard (least productive)"))) |>
    asF()
fig7 <- ggplot(sw, aes(yield, rel, colour = regime)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    geom_vline(xintercept = r$target, colour = MUTED, linetype = "22", linewidth = 0.4) +
    facet_wrap(~metric) +
    regime_scale +
    scale_x_log10(labels = label_log()) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = "In this ecosystem, BH_P does not buy any protection",
         subtitle = paste0("Each rule swept over fishing intensity; the dashed line is the calibrated yield used in the other figures.\n",
             "BH_P fishes Gurnard at 1e-5 per year - effectively not at all - yet Gurnard ends up slightly worse off than under a\n",
             "fixed F (right). What reaches it is not the hook but the food web. And the species that limits the fishery here is\n",
             "Cod, the most productive one, which BH_P fishes hardest - so BH_P is also the worst of the three on the usual\n",
             "most-depleted test (left)."),
         x = expression(total~yield~(g~yr^-1)), y = "biomass, % of unfished",
         caption = paste0("This is not evidence against Law & Plank. NS_params is not their setting: twelve strongly coupled, well-studied stocks, none rare\n",
             "enough to be at risk of collapse. Their result concerns assembled ecosystems containing genuinely rare species.")) +
    theme_bh()

## Cover page ---------------------------------------------------------------
sig <- function(x, d = 3) formatC(x, format = "g", digits = d)
# each row: style, text, and (for key/value rows) the value column
rows <- list(
  c("t", "Balanced harvesting in mizer", ""),
  c("s", "Size-independent fishing mortality, after Law & Plank (2023)", ""),
  c("rule", "", ""),
  c("h", "The three rules", ""),
  c("m", "fixed        F_i(t)  =  F_i", ""),
  c("m", "BH_P         F_i(t)  =  c_P P_i(t)", ""),
  c("m", "BH_P/B       F_i(t)  =  c_P/B P_i(t) / B_i(t)", ""),
  c("gap", "", ""),
  c("b", "In all three, fishing mortality is independent of body size. Every species enters one mixed fishery at a", ""),
  c("b", "single entry mass w_f, and all fish above it are caught at the same rate. The balancing is across species,", ""),
  c("b", "not across sizes. P_i and B_i are the somatic production rate and the biomass over the harvested range.", ""),
  c("gap", "", ""),
  c("h", "This run", ""),
  c("k", "Model", paste0("North Sea (mizer NS_params ", packageVersion("mizer"), "), 12 species")),
  c("k", "Entry mass", paste0("w_f = ", r$w_f, " g, the value used in the paper")),
  c("k", "Starting state", paste0("unfished steady state, then ", r$t_max, " years of fishing")),
  c("k", "Implementation", "z_i ~ U(0.5, 1.5) per species, fixed in time, as in their Fig. 5"),
  c("k", "Calibration", paste0("all three tuned to the same year-", r$t_max, " total yield, ", sig(r$target), " g/yr")),
  c("km", "", paste0("F = 0.2 z_i     c_P = ", sig(r$const[["BHP"]]), "     c_P/B = ", sig(r$const[["BHPB"]]))),
  c("k", "Caution", "c_P has dimensions area/mass, so its value belongs to this model's biomass"),
  c("k", "", "units and does not carry over to another."),
  c("gap", "", ""),
  c("h", "Figures", ""),
  c("n", "1", "Fishing mortality does not depend on body size"),
  c("n", "2", "Balancing is across species, not across sizes"),
  c("n", "3", "Under BH_P the fishing mortality is a feedback, not a setting"),
  c("n", "4", "Stock response over 50 years of fishing, at equal total yield"),
  c("n", "5", "Why BH_P takes so little from the unproductive species"),
  c("n", "6", "Fishing the big fish lifts the middle of the spectrum"),
  c("n", "7", "In this ecosystem, BH_P does not buy any protection"),
  c("gap", "", ""),
  c("c", "Law, R. & Plank, M. J. (2023). Fishing for biodiversity by balanced harvesting. Fish and Fisheries.", ""),
  c("c", "doi:10.1111/faf.12705    -    Generated by bh_run.R and bh_figures.R in this project.", ""))

# style -> c(height, size_mm, x_label, x_value, bold, mono, colour)
st <- list(
  t    = list(h = 2.4, sz = 6.4, x = 0,   xv = NA, face = 2, fam = "sans", col = INK),
  s    = list(h = 2.0, sz = 3.9, x = 0,   xv = NA, face = 1, fam = "sans", col = INK2),
  h    = list(h = 1.9, sz = 3.8, x = 0,   xv = NA, face = 2, fam = "sans", col = INK),
  b    = list(h = 1.1, sz = 3.1, x = 0,   xv = NA, face = 1, fam = "sans", col = INK2),
  m    = list(h = 1.2, sz = 3.1, x = 1.5, xv = NA, face = 1, fam = "mono", col = INK),
  k    = list(h = 1.1, sz = 3.1, x = 0,   xv = 13, face = 1, fam = "sans", col = INK2),
  km   = list(h = 1.3, sz = 3.1, x = 0,   xv = 13, face = 1, fam = "mono", col = INK),
  n    = list(h = 1.1, sz = 3.1, x = 0.3, xv = 2,  face = 1, fam = "sans", col = INK2),
  c    = list(h = 1.0, sz = 2.7, x = 0,   xv = NA, face = 1, fam = "sans", col = MUTED),
  gap  = list(h = 0.8, sz = 0.1, x = 0,   xv = NA, face = 1, fam = "sans", col = SURFACE),
  rule = list(h = 1.4, sz = 0.1, x = 0,   xv = NA, face = 1, fam = "sans", col = SURFACE))

hts <- vapply(rows, function(x) st[[x[1]]]$h, 1)
y   <- -(cumsum(hts) - hts / 2)
d   <- data.frame(y = y, style = vapply(rows, `[`, "", 1),
                  lab = vapply(rows, `[`, "", 2), val = vapply(rows, `[`, "", 3))
get1 <- function(s, f) vapply(s, function(k) st[[k]][[f]], if (f %in% c("fam","col")) "" else 1)

cover <- ggplot(d) +
    annotate("segment", x = 0, xend = 66, colour = "#d8d7d1", linewidth = 0.4,
             y = y[d$style == "rule"], yend = y[d$style == "rule"]) +
    geom_text(aes(x = get1(style, "x"), y = y, label = lab, size = get1(style, "sz"),
                  fontface = get1(style, "face"), family = get1(style, "fam"),
                  colour = get1(style, "col")), hjust = 0) +
    geom_text(data = subset(d, val != ""),
              aes(x = get1(style, "xv"), y = y, label = val, size = get1(style, "sz"),
                  family = get1(style, "fam"), colour = get1(style, "col")), hjust = 0) +
    scale_size_identity() + scale_colour_identity() +
    scale_x_continuous(limits = c(0, 66)) +
    scale_y_continuous(limits = c(min(y) - 1, 0)) +
    theme_void() +
    theme(plot.background = element_rect(fill = SURFACE, colour = NA),
          legend.position = "none", plot.margin = margin(26, 26, 22, 28))

## Write the PDF ------------------------------------------------------------
pdf("balanced_harvesting.pdf", width = 9.6, height = 6.4, onefile = TRUE)
print(cover); for (i in 1:7) print(get(paste0("fig", i)))
invisible(dev.off())
message("wrote balanced_harvesting.pdf")
