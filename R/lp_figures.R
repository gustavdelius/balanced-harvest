# ---------------------------------------------------------------------------
# Figures 2-6 of Law & Plank (2023).
#
# Figure 2a is a scatter of biomass against production rate from the Ecopath
# model of the West Scotland shelf (Alexander et al. 2015).  That is empirical
# data, not model output, and is not reproduced here; Fig. 2b, c are.
# ---------------------------------------------------------------------------

source("R/lp_harvest.R")
library(ggplot2)
library(patchwork)

RULES  <- c("fixed", "BHP", "BHPB")
RULE_LAB <- c(fixed = "F[i]*': constant'",
              BHP   = "F[i]*': '*BH[P]",
              BHPB  = "F[i]*': '*BH[P/B]")

theme_lp <- function(base = 9) {
    theme_bw(base_size = base) +
        theme(panel.grid = element_blank(),
              strip.background = element_blank(),
              strip.text = element_text(size = base, hjust = 0.5),
              plot.tag = element_text(size = base, face = "plain"),
              plot.tag.position = c(0.03, 0.97),
              legend.position = "none")
}

# scales::breaks_log() returns fractional powers of ten when the data span
# little more than a decade, which label_log() then renders as "10^0.477".
# Always break on whole decades instead.
decade_breaks <- function(n = 5) function(x) {
    x <- x[is.finite(x) & x > 0]
    if (!length(x)) return(numeric(0))
    lo <- floor(log10(min(x)));  hi <- ceiling(log10(max(x)))
    10^seq(lo, hi, by = max(1, ceiling((hi - lo) / n)))
}

scale_log <- function(axis = "y", ...) {
    f <- if (axis == "y") scale_y_log10 else scale_x_log10
    f(breaks = decade_breaks(), labels = scales::label_log(), ...)
}

## --- Growth trajectories, Eq. (C.1) ----------------------------------------
# dw/dt = eps_i(w) g_i(w) = mizer's e_growth, integrated up from egg size.
lp_growth_curve <- function(params, n = initialN(params),
                            n_pp = initialNResource(params)) {
    e <- getEGrowth(params, n = n, n_pp = n_pp)
    w <- w(params)
    w_max <- species_params(params)$w_max
    do.call(rbind, lapply(seq_len(nrow(e)), function(i) {
        age <- cumsum(c(0, diff(w)) / pmax(e[i, ], .Machine$double.xmin))
        keep <- w <= w_max[i] & is.finite(age)
        data.frame(species = i, age = age[keep], w = w[keep])
    }))
}

## --- Figure 2 --------------------------------------------------------------
# (b) biomass against production rate, (c) biomass against maximum body mass.
# Both measured over the whole size range of each species, as in the paper.
lp_fig2 <- function(params) {
    n  <- initialN(params)
    eg <- getEGrowth(params)
    d <- data.frame(
        sp    = seq_len(nrow(n)),
        w_max = species_params(params)$w_max,
        B     = lp_biomass(params, n, 1),
        P     = lp_production(params, n, eg, 1))

    b <- ggplot(d, aes(P, B)) +
        geom_text(aes(label = sp), size = 2.6) +
        scale_log("x") + scale_log("y") +
        labs(x = expression("production rate (g m"^-2*" yr"^-1*")"),
             y = expression("biomass (g m"^-2*")"), tag = "(b)") +
        theme_lp()
    cc <- ggplot(d, aes(w_max, B)) +
        geom_text(aes(label = sp), size = 2.6) +
        scale_log("x") + scale_log("y") +
        labs(x = "max body mass (g)",
             y = expression("biomass (g m"^-2*")"), tag = "(c)") +
        theme_lp()
    list(plot = b + cc, data = d)
}

## --- Figure 3 --------------------------------------------------------------
# Rows: (a-c) biomass time series, (d-f) F_i in year 0 and year 50,
# (g-i) ratio of year-50 to year-0 biomass.  Columns: the three rules.
lp_fig3 <- function(tracks, t_max = LP_FISHING$t_max) {
    tr <- do.call(rbind, lapply(RULES, function(r)
        transform(tracks[[r]], rule = r)))
    tr$rule <- factor(tr$rule, RULES, RULE_LAB[RULES])
    tr$sp   <- as.integer(tr$species)

    ends <- subset(tr, time %in% c(0, t_max))
    rel <- merge(subset(tr, time == t_max, c(rule, sp, B_total)),
                 subset(tr, time == 0,     c(rule, sp, B_total)),
                 by = c("rule", "sp"), suffixes = c("_end", "_start"))
    rel$ratio <- rel$B_total_end / rel$B_total_start

    p1 <- ggplot(tr, aes(time, B_total, group = sp)) +
        geom_line(linewidth = 0.25) +
        facet_wrap(~rule, nrow = 1, labeller = label_parsed) +
        scale_log("y") +
        labs(x = "time (yr)", y = expression("total biomass (g m"^-2*")")) +
        theme_lp()
    p2 <- ggplot(subset(ends, F > 0), aes(sp, F, shape = factor(time))) +
        geom_point(size = 1.4, fill = "white") +
        facet_wrap(~rule, nrow = 1, labeller = label_parsed) +
        scale_shape_manual(values = c(21, 19)) +
        scale_log("y") +
        labs(x = "species", y = expression(F[i]*" (yr"^-1*")")) +
        theme_lp()
    p3 <- ggplot(rel, aes(sp, ratio)) +
        geom_hline(yintercept = 1, linetype = 3, linewidth = 0.3) +
        geom_point(size = 1.2) +
        facet_wrap(~rule, nrow = 1, labeller = label_parsed) +
        scale_log("y") +
        labs(x = "species", y = "relative biomass") +
        theme_lp()
    p1 / p2 / p3
}

## --- Figure 4 --------------------------------------------------------------
# Yield against production rate over the harvested size range, at year 50.
# Dotted lines: constant exploitation ratio E = c_P/B.  Dashed line in the
# BH_P panel: slope 1 + alpha, where B ~ P^alpha across the assemblage.
lp_fig4 <- function(tracks, params, t_max = LP_FISHING$t_max,
                    E = LP_FISHING$c_PB, track_species = NULL) {
    end <- do.call(rbind, lapply(RULES, function(r)
        transform(subset(tracks[[r]], time == t_max), rule = r)))
    end$rule <- factor(end$rule, RULES, RULE_LAB[RULES])
    end$sp <- as.integer(end$species)
    end <- subset(end, Y > 0 & P > 0)

    # alpha from B ~ P^alpha over the harvested range in the unfished state
    params <- lp_set_fishing(params, "none")   # defines the harvested range
    n0 <- initialN(params)
    d0 <- data.frame(B = lp_biomass(params, n0),
                     P = lp_production(params, n0, getEGrowth(params)))
    d0 <- subset(d0, B > 0 & P > 0)
    alpha <- unname(coef(lm(log(B) ~ log(P), d0))[2])

    rng <- range(end$P)
    # dotted isocline of constant exploitation ratio E, in the fixed and
    # BH_P/B panels only; dashed line of slope 1 + alpha in the BH_P panel.
    iso <- do.call(rbind, lapply(c("fixed", "BHPB"), function(r)
        data.frame(P = rng, Y = E * rng,
                   rule = factor(RULE_LAB[[r]], levels(end$rule)),
                   row.names = NULL)))
    bh <- subset(end, rule == RULE_LAB["BHP"])
    k  <- exp(mean(log(bh$Y) - (1 + alpha) * log(bh$P)))
    slp <- data.frame(P = rng, Y = k * rng^(1 + alpha),
                      rule = factor(RULE_LAB[["BHP"]], levels(end$rule)),
                      row.names = NULL)

    traj <- NULL
    if (!is.null(track_species)) {
        traj <- do.call(rbind, lapply(RULES, function(r)
            transform(subset(tracks[[r]], species == track_species & P > 0 & Y > 0),
                      rule = r)))
        traj$rule <- factor(traj$rule, RULES, RULE_LAB[RULES])
    }

    p <- ggplot(end, aes(P, Y)) +
        geom_line(data = iso, aes(P, Y), linetype = 3, linewidth = 0.3,
                  inherit.aes = FALSE) +
        geom_line(data = slp, aes(P, Y), linetype = 2, linewidth = 0.3,
                  inherit.aes = FALSE)
    if (!is.null(traj))
        p <- p + geom_path(data = traj, aes(P, Y), linewidth = 0.4,
                           colour = "grey45",
                           arrow = arrow(length = unit(0.06, "inches"),
                                         type = "closed")) +
            geom_point(data = subset(traj, time == 0), aes(P, Y),
                       shape = 21, fill = "white", colour = "grey45",
                       size = 1.6)
    p + geom_text(aes(label = sp), size = 2.6) +
        facet_wrap(~rule, nrow = 1, labeller = label_parsed) +
        scale_log("x") + scale_log("y") +
        # the reference lines run the full width of the P axis, which would
        # otherwise drag the shared y axis far below the data
        coord_cartesian(ylim = range(c(end$Y, traj$Y))) +
        labs(x = expression("production rate (g m"^-2*" yr"^-1*")"),
             y = expression("yield (g m"^-2*" yr"^-1*")")) +
        theme_lp() -> p
    list(plot = p, alpha = alpha)
}

## --- Figure 5 --------------------------------------------------------------
# Growth trajectories in three independently assembled ecosystems.
lp_fig5 <- function(params_list) {
    d <- do.call(rbind, lapply(seq_along(params_list), function(k)
        transform(lp_growth_curve(params_list[[k]]),
                  assemblage = paste("Assemblage", k))))
    ggplot(d, aes(age, w, group = species)) +
        geom_line(linewidth = 0.25) +
        facet_wrap(~assemblage, nrow = 1) +
        scale_log("y") +
        coord_cartesian(xlim = c(0, 25)) +
        labs(x = "age (yr)", y = "body mass (g)") +
        theme_lp()
}

## --- Figure 6 --------------------------------------------------------------
# Biomass time series for three assemblages under the three rules.
lp_fig6 <- function(tracks_list) {
    d <- do.call(rbind, lapply(seq_along(tracks_list), function(k)
        do.call(rbind, lapply(RULES, function(r)
            transform(tracks_list[[k]][[r]], rule = r,
                      assemblage = paste("Assemblage", k))))))
    d$rule <- factor(d$rule, RULES, RULE_LAB[RULES])
    d$sp <- as.integer(d$species)
    ggplot(d, aes(time, B_total, group = sp)) +
        geom_line(linewidth = 0.25) +
        facet_grid(assemblage ~ rule, labeller = labeller(rule = label_parsed)) +
        scale_log("y") +
        labs(x = "time (yr)", y = expression("total biomass (g m"^-2*")")) +
        theme_lp()
}
