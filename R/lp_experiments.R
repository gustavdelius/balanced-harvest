# ---------------------------------------------------------------------------
# Machinery for the robustness investigation of docs/robustness.md.
#
# The central idea (robustness.md section 5.1) is to stop comparing rules at a
# single calibrated point and instead sweep each rule's intensity constant,
# giving a yield-versus-biodiversity frontier per rule.  A rule whose frontier
# dominates does so whatever criterion is used to match yields, which removes
# the arbitrary "equal yield after 50 years" choice entirely.
# ---------------------------------------------------------------------------

source("R/lp_harvest.R")

## --- Choosing a comparable intensity for each rule -------------------------
# Every rule is F_i = c z'_i g_i.  Setting c so that the *initial* total yield
# matches that of fixed fishing at F_ref puts all rules on the same footing
# before any sweeping, and costs nothing: Y(0) = c sum_i z_i g_i B_i(0).
lp_base_const <- function(params, rule, zf = NULL, theta = NULL, alloc = NULL,
                          F_ref = LP_FISHING$F_base, w_f = LP_FISHING$w_f,
                          w_measure = NULL) {
    if (is.null(zf)) zf <- rep(1, nrow(species_params(params)))
    p <- lp_set_fishing(params, rule, 1, zf, w_f = w_f,
                        theta = theta, alloc = alloc, w_measure = w_measure)
    n <- initialN(p)
    B <- lp_biomass(p, n)
    g <- lp_alloc(p, n, getEGrowth(p))
    g[!is.finite(g)] <- 0
    denom <- sum(zf * g * B)
    if (denom <= 0) return(NA_real_)
    F_ref * sum(B) / denom
}

# The year-0 allocation of another rule, for the `frozen` controls.
lp_frozen_alloc <- function(params, rule, w_f = LP_FISHING$w_f) {
    p <- lp_set_fishing(params, rule, 1, w_f = w_f)
    g <- lp_alloc(p, initialN(p), getEGrowth(p))
    g[!is.finite(g)] <- 0
    g
}

# The body mass above which a given fraction of community biomass sits.  Used
# to give variants with very different size structures a comparable fishery:
# matching w_f on body mass alone is meaningless when contesting the plankton
# crushes juvenile growth and leaves nothing near w_max (docs/robustness.md
# section 3.6).
lp_wf_for_fraction <- function(params, frac) {
    n <- initialN(params); w <- w(params); dw <- dw(params)
    bw <- colSums(sweep(n, 2, w * dw, "*"))
    cum <- rev(cumsum(rev(bw))) / sum(bw)
    w[which.min(abs(cum - frac))]
}

lp_fraction_above <- function(params, w_f) {
    n <- initialN(params); w <- w(params); dw <- dw(params)
    bw <- colSums(sweep(n, 2, w * dw, "*"))
    sum(bw[w >= w_f]) / sum(bw)
}

## --- Outcome measures ------------------------------------------------------
# All biodiversity measures are relative to the unfished control at the same
# time, not to year 0, because the assembled state is only a quasi-equilibrium
# and drifts on its own (robustness.md section 3.3).
shannon <- function(B) {
    p <- B[B > 0] / sum(B[B > 0])
    -sum(p * log(p))
}

lp_metrics <- function(track, control, t_max = LP_FISHING$t_max) {
    e  <- subset(track,   time == t_max)
    s0 <- subset(track,   time == 0)
    c50 <- subset(control, time == t_max)
    e   <- e[match(s0$species, e$species), ]
    c50 <- c50[match(s0$species, c50$species), ]

    rel  <- e$B_total / s0$B_total                 # as the paper reports it
    relc <- e$B_total / c50$B_total                # against no fishing
    relc[!is.finite(relc)] <- 0

    # cumulative yield over the 50 years, by trapezoid on the saved times
    ty <- aggregate(Y ~ time, track, sum)
    ty <- ty[order(ty$time), ]
    cum <- sum(diff(ty$time) * (head(ty$Y, -1) + tail(ty$Y, -1)) / 2)

    data.frame(
        yield       = sum(e$Y),
        yield_cum   = cum,
        min_rel     = min(rel),
        min_rel_ctrl = min(relc),
        n_below_10pct_ctrl = sum(relc < 0.1),
        n_below_1pct_ctrl  = sum(relc < 0.01),
        # overall distortion of the assemblage, in log space
        rms_log_ctrl = sqrt(mean(log(pmax(relc, 1e-12))^2)),
        shannon_ctrl = shannon(e$B_total) - shannon(c50$B_total),
        row.names = NULL)
}

## --- One rule, swept over intensity ----------------------------------------
LP_MULTS <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)

lp_frontier <- function(params, rule, control, zf = NULL, theta = NULL,
                        alloc = NULL, mults = LP_MULTS,
                        F_ref = LP_FISHING$F_base, w_f = LP_FISHING$w_f,
                        w_measure = NULL,
                        dt = 0.01, t_save = 5, label = rule) {
    if (is.null(zf)) zf <- rep(1, nrow(species_params(params)))
    base <- lp_base_const(params, rule, zf, theta, alloc, F_ref, w_f, w_measure)
    do.call(rbind, lapply(mults, function(m) {
        sim <- lp_harvest(params, rule, base * m, zf, dt = dt,
                          t_save = t_save, theta = theta, alloc = alloc,
                          w_f = w_f, w_measure = w_measure)
        cbind(data.frame(rule = label, mult = m, const = base * m),
              lp_metrics(lp_track(sim), control))
    }))
}

# The unfished control must be run at the same dt as the sweep it is compared
# against, so that numerical differences cancel.
lp_control <- function(params, zf = NULL, dt = 0.01, t_save = 5,
                       w_f = LP_FISHING$w_f) {
    if (is.null(zf)) zf <- rep(1, nrow(species_params(params)))
    lp_track(lp_harvest(params, "none", 0, zf, dt = dt, t_save = t_save,
                        w_f = w_f))
}

## --- Reading a frontier at a common yield ----------------------------------
# Log-log interpolation of a metric along a rule's frontier, so that rules can
# be compared at a yield none of them was run at exactly.
lp_at_yield <- function(front, metric, y) {
    f <- front[order(front$mult), ]
    ok <- f$yield > 0 & is.finite(f[[metric]])
    f <- f[ok, ]
    if (nrow(f) < 2) return(NA_real_)
    # Yield is not always monotone in intensity: rules without a stabilising
    # feedback overfish into declining yield, so the frontier can double back.
    # Interpolating across the turnover would mix the two branches, so keep
    # only the ascending branch up to the yield maximum.
    f <- f[seq_len(which.max(f$yield)), , drop = FALSE]
    if (nrow(f) < 2) return(NA_real_)
    # NA rather than an extrapolation when the rule cannot reach this yield at
    # any intensity - that is a result, not missing data.
    if (y < min(f$yield) || y > max(f$yield)) return(NA_real_)
    v <- f[[metric]]
    if (all(v > 0)) exp(approx(log(f$yield), log(v), log(y))$y)
    else approx(log(f$yield), v, log(y))$y
}

# The largest terminal yield a rule reaches anywhere on its frontier.
lp_max_yield <- function(front) max(front$yield, na.rm = TRUE)

## --- Phase 4: the rule as it would actually be implemented -----------------
# The rules above recompute F_i at every time step from perfectly known P_i and
# B_i.  A real harvest control rule is updated periodically from survey
# estimates.  lp_harvest_periodic() projects in segments of `interval` years:
# at the start of each segment it observes the state, optionally corrupts the
# observation with log-normal error of log-scale s.d. `sigma` (mean-preserving),
# freezes F_i at that value, and projects the segment.
#
#   g_hat_i = g_i(t_k) * exp(eta_ik),   eta_ik ~ N(-sigma^2/2, sigma^2)
#
# `interval = 0` means continuous updating, i.e. the idealised rule used above.
lp_harvest_periodic <- function(params, rule, const, zf = NULL, interval = 0,
                                sigma = 0, t_max = LP_FISHING$t_max,
                                dt = 0.01, t_save = 5, seed = NULL) {
    if (is.null(zf)) zf <- rep(1, nrow(species_params(params)))
    if (interval <= 0 && sigma == 0)
        return(lp_track(lp_harvest(params, rule, const, zf, t_max = t_max,
                                   dt = dt, t_save = t_save)))
    if (!is.null(seed)) set.seed(seed)
    p <- params
    t0 <- 0
    out <- list()
    while (t0 < t_max - 1e-9) {
        len <- min(interval, t_max - t0)
        pp <- lp_set_fishing(p, rule, const, zf)
        g  <- lp_alloc(pp, initialN(pp), getEGrowth(pp))
        g[!is.finite(g)] <- 0
        if (sigma > 0)
            g <- g * exp(rnorm(length(g), -sigma^2 / 2, sigma))
        pf  <- lp_set_fishing(p, "frozen", const, zf, alloc = g)
        sim <- project(pf, t_max = len, dt = dt, t_save = min(t_save, len),
                       progress_bar = FALSE)
        tr <- lp_track(sim)
        tr$time <- tr$time + t0
        if (t0 > 0) tr <- subset(tr, time > t0 + 1e-9)
        out[[length(out) + 1]] <- tr
        p  <- finalParams(sim)
        t0 <- t0 + len
    }
    do.call(rbind, out)
}
