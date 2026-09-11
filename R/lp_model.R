# ---------------------------------------------------------------------------
# The Law & Plank (2023) ecosystem model, expressed as a mizer MizerParams
# object.  Each block below names the equation of the paper it implements.
#
# The mapping is exact except where noted in `docs/mapping.md`.  The main
# points are:
#   * no satiation and no metabolic cost, so E = K * encounter  (Eq. A.1a)
#   * intrinsic mortality is larval + background, and the background part is
#     food-dependent, hence a custom `Mort` rate function        (Eqs A.5-A.7)
#   * the plankton follow a logistic with immigration rather than mizer's
#     semichemostat, hence a custom resource dynamics            (Eq. A.11)
#   * reproduction has no stock-recruitment relationship         (Eq. A.9)
# ---------------------------------------------------------------------------

library(mizer)
source("R/lp_constants.R")

## --- Species parameters ----------------------------------------------------
# A life-history template (Fig. 1) in which species differ only in maximum
# body mass, larval death rate, and the optional activity factor z_i.
lp_species_params <- function(w_max, mu_egg, z = 1,
                              species = NULL, fish = LP_FISH) {
    n <- length(w_max)
    if (is.null(species)) species <- as.character(seq_len(n))
    data.frame(
        species = species,
        w_max   = w_max,
        w_min   = fish$w_egg,
        # w_mat = w_max/10, and w_mat25 chosen so that mizer's maturity
        # exponent U = log(3)/log(w_mat/w_mat25) equals rho_m      Eq. (A.8a)
        w_mat   = w_max * fish$w_mat_rat,
        w_mat25 = w_max * fish$w_mat_rat / 3^(1 / fish$rho_m),
        # repro_prop = (w/w_max)^(m - n) must equal (w/w_max)^rho_inf Eq. (A.8b)
        n = 2/3, m = 2/3 + fish$rho_inf,
        # box feeding kernel over [beta - 3 sigma, beta + 3 sigma]  Eq. (A.4)
        pred_kernel_type = "box",
        ppmr_min = exp(fish$beta - 3 * fish$sigma),
        ppmr_max = exp(fish$beta + 3 * fish$sigma),
        # search volume A_i w^alpha; the 1/(6 sigma) normalisation of the box
        # kernel is folded into gamma because mizer's box kernel has height 1
        gamma = z * fish$A / (6 * fish$sigma),
        q     = fish$alpha_q,
        # no satiation (h = Inf => feeding level 0) and no metabolic cost
        h = Inf, ks = 0, k = 0,
        alpha = fish$K,                       # conversion efficiency K
        interaction_resource = 1,             # theta_{i0} = 1
        # mizer's RDI carries a factor 1/2 for the sex ratio, which the paper
        # does not, so erepro = 2 * eps_R                          Eq. (A.9)
        erepro = 2 * fish$eps_R,
        # carried through for the mortality functions
        mu_egg = mu_egg,
        z      = z,
        stringsAsFactors = FALSE
    )
}

# theta: cannibalism stronger than between-species predation    (Appendix C)
lp_interaction <- function(n, fish = LP_FISH) {
    m <- matrix(fish$theta_ij, nrow = n, ncol = n)
    diag(m) <- fish$theta_ii
    m
}

## --- Intrinsic mortality ---------------------------------------------------
# Larval component: a reverse sigmoid, large on eggs, vanishing by w_L.
# Static in time, so it goes into mizer's ext_mort slot.                (A.6)
lp_larval_mort <- function(params, fish = LP_FISH) {
    sp <- species_params(params)
    outer(sp$z * sp$mu_egg, w(params),
          function(mu, w) mu / (1 + (w / fish$w_L)^fish$rho_L))
}

# Background component.  Eq. (A.7) as printed is
#
#   mu_b,i(w,t) = mu_b0 (w/w_0)^(-xi) * g_i(w,t)/g_i(w_0,t)            "product"
#
# but the prose describing it says mu_b is "proportional to the mass-specific
# needs for metabolism, RELATIVE TO the mass-specific rate at which food
# becomes available", i.e.
#
#   mu_b,i(w,t) = mu_b0 (w/w_0)^(-xi) / [g_i(w,t)/g_i(w_0,t)]          "ratio"
#
# The two differ drastically: see docs/mu_b.md.  `mu_b_form` selects between
# them.  Depends on the mass-specific intake rate g_i(w, t) either way, so it
# changes as the model runs.  The mass-specific intake g_i(w,t) is
# K * encounter_i(w,t) / w, and the factor K cancels in the ratio.
#
# The species- and size-dependent prefactor z_i mu_b0 (w/w_0)^(-xi) is static,
# so lp_params() precomputes it into other_params(params)$lp_mu_b_pre.
lp_background_mort <- function(params, encounter) {
    pre <- other_params(params)$lp_mu_b_pre
    w   <- params@w
    # g_i(w, t) / g_i(w_0, t) with g mass-specific; the factor K cancels
    g_rel <- sweep(sweep(encounter, 2, w, "/"), 1, encounter[, 1] / w[1], "/")
    if (identical(other_params(params)$lp_mu_b_form, "ratio")) {
        # mu_b proportional to metabolic need DIVIDED BY food availability
        g_rel <- 1 / g_rel
    }
    g_rel[!is.finite(g_rel)] <- 0
    pre * g_rel
}

lp_mu_b_prefactor <- function(params, fish = LP_FISH) {
    sp <- species_params(params)
    outer(sp$z * fish$mu_b0, params@w / params@w[1], function(a, r) a * r^(-fish$xi))
}

# Total mortality: predation + fishing + the larval term (which lives in
# mizer's ext_mort slot, `params@mu_b` - not to be confused with the paper's
# mu_b, which is the background term added here).
#
# mizer's plumbing does not hand the encounter rate to the Mort function, so it
# has to be recomputed, one extra convolution per time step.  Caching it
# between the Encounter and Mort calls would be unsafe because getMort() can be
# called on its own.
lpMort <- function(params, n, n_pp, n_other, t = 0, f_mort, pred_mort, ...) {
    encounter <- mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other,
                                t = t)
    pred_mort + params@mu_b + f_mort + lp_background_mort(params, encounter)
}

## --- Feeding level ---------------------------------------------------------
# The paper has no maximum intake rate, so the feeding level is identically
# zero.  h = Inf would give this too, but an explicit zero avoids Inf/Inf.
lpFeedingLevel <- function(params, n, n_pp, n_other, t = 0, encounter, ...) {
    encounter * 0
}

## --- Diffusion -------------------------------------------------------------
# Term (f) of Eq. (A.1) carries epsilon_i(x) G_i(x, t), i.e. only the part of
# the food that goes to somatic growth diffuses.  mizer's built-in diffusion
# omits that factor, so multiply it back in.  epsilon_i = 1 - psi.
lpDiffusion <- function(params, n, n_pp, n_other, t = 0, feeding_level, ...) {
    (1 - params@psi) *
        mizerDiffusion(params, n, n_pp, n_other, t = t,
                       feeding_level = feeding_level, ...)
}

## --- Plankton dynamics -----------------------------------------------------
# Logistic growth, predation by fish, and immigration:               Eq. (A.11)
#
#   dn/dt = I + (r - d) n - (r/a) n^2
#
# The paper integrated this with explicit Euler at dt = 0.002 yr.  That is a
# Riccati equation with coefficients that are constant over a time step, so we
# use its closed-form solution instead: same equation, same fixed points, but
# unconditionally stable, which frees the time step from the stiffness of the
# fastest-dividing cells (r reaches ~320/yr at 1e-10 g).  See
# `tests/test_convergence.R` for the check that this agrees with Euler at
# dt = 0.002.
#
# Writing n' = k (n+ - n)(n - n-) with n+ > 0 >= n- the roots of
# k n^2 - (r - d) n - I = 0, the exact step is
#
#   v = (n+ - n)/(n - n-) * exp(-k (n+ - n-) dt),   n_new = (n+ + n- v)/(1 + v)
#
lpResource <- function(params, n, n_pp, n_other, rates, t, dt,
                       resource_rate, resource_capacity, ...) {
    idx <- other_params(params)$lp_resource_idx
    imm <- other_params(params)$lp_immigration[idx]
    r   <- resource_rate[idx]
    a   <- resource_capacity[idx]
    u   <- n_pp[idx]

    k    <- r / a
    b    <- r - rates$resource_mort[idx]
    disc <- sqrt(b * b + 4 * k * imm)
    n_hi <- (b + disc) / (2 * k)
    n_lo <- (b - disc) / (2 * k)

    v <- (n_hi - u) / (u - n_lo) * exp(-disc * dt)
    out <- n_pp
    out[idx] <- (n_hi + n_lo * v) / (1 + v)
    out[out < 0] <- 0
    out
}

# The explicit Euler step the paper used, kept for the convergence check.
lpResourceEuler <- function(params, n, n_pp, n_other, rates, t, dt,
                            resource_rate, resource_capacity, ...) {
    idx <- other_params(params)$lp_resource_idx
    imm <- other_params(params)$lp_immigration[idx]
    r <- resource_rate[idx]
    a <- resource_capacity[idx]
    u <- n_pp[idx]
    out <- n_pp
    out[idx] <- u + dt * (r * u * (1 - u / a) -
                              rates$resource_mort[idx] * u + imm)
    out[out < 0] <- 0
    out
}

## --- Assembling the MizerParams object -------------------------------------
lp_params <- function(w_max, mu_egg, z = rep(1, length(w_max)),
                      species = NULL,
                      fish = LP_FISH, plankton = LP_PLANKTON,
                      numerics = LP_NUMERICS,
                      grid_w_max = LP_ASSEMBLY$w_max_range[2],
                      mu_b_form = LP_FISH$mu_b_form) {
    stopifnot(all(w_max <= grid_w_max))
    sp <- lp_species_params(w_max, mu_egg, z, species, fish)

    # Grid: a fixed log step dx over the whole possible range, so that the
    # grid does not move as species are added or removed during assembly.
    no_w <- round(log(grid_w_max / fish$w_egg) / numerics$dx) + 1

    params <- newMultispeciesParams(
        sp,
        interaction = lp_interaction(nrow(sp), fish),
        no_w      = no_w,
        min_w     = fish$w_egg,
        max_w     = grid_w_max,
        min_w_pp  = plankton$w_min,
        n         = 2/3,
        lambda    = plankton$lambda,
        w_pp_cutoff = plankton$w_max,
        resource_dynamics = "lpResource",
        RDD       = "noRDD",
        info_level = 0
    )

    # Predation-driven diffusion, term (f) of Eq. (A.1)
    use_predation_diffusion(params) <- TRUE

    # Plankton rates, converted from the paper's log-densities u(x) to mizer's
    # per-mass densities phi(w) = u(x)/w.                     Eqs (A.12)-(A.14)
    wf  <- w_full(params)
    on  <- wf <= plankton$w_max
    cap <- ifelse(on, plankton$a0 * (wf / plankton$w_a)^(1 - plankton$lambda) / wf, 0)
    imm <- ifelse(on, plankton$I0 * (wf / plankton$w_I)^(1 - plankton$lambda) / wf, 0)
    rate <- ifelse(on, plankton$r0 * wf^(-plankton$rho), 0)
    resource_capacity(params) <- cap
    resource_rate(params)     <- rate
    initialNResource(params)  <- cap

    other_params(params)$lp_immigration   <- imm
    other_params(params)$lp_resource_idx  <- which(on)
    other_params(params)$lp_fish          <- fish
    other_params(params)$lp_plankton      <- plankton

    # Larval mortality (Eq. A.6); the background part is added by lpMort.
    ext_mort(params) <- lp_larval_mort(params, fish)
    other_params(params)$lp_mu_b_pre  <- lp_mu_b_prefactor(params, fish)
    other_params(params)$lp_mu_b_form <- mu_b_form

    params <- setRateFunction(params, "FeedingLevel", "lpFeedingLevel")
    params <- setRateFunction(params, "Mort",         "lpMort")
    params <- setRateFunction(params, "Diffusion",    "lpDiffusion")
    params
}

## --- Initial abundance of an invading species ------------------------------
# "The size spectrum was started at a low density on a power law, the egg
# density being 0.002 m^-2" (Appendix B).  The quoted density is the
# log-density u(x_0), so phi(w_egg) = 0.002 / w_egg.
lp_invader_n <- function(params, idx, assembly = LP_ASSEMBLY) {
    w <- w(params)
    n0 <- assembly$n_egg_init / w[1]
    n <- n0 * (w / w[1])^(-assembly$invader_slope)
    n[w > species_params(params)$w_max[idx]] <- 0
    n
}
