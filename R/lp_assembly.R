# ---------------------------------------------------------------------------
# Appendix B: sequential assembly of the model ecosystems.
#
# Start from plankton alone and add one fish species at a time, letting the
# system relax for 50 years after each invasion and removing species whose
# biomass has fallen below an extinction threshold.  Stop at 15 species or
# 40 attempts, whichever comes first.
#
# The paper does not report the random seed or the realised life histories of
# its ecosystems, so this reproduces the procedure, not the particular
# assemblage.  Every assemblage that satisfies the four selection criteria at
# the end of Appendix B is an equally good starting point.
# ---------------------------------------------------------------------------

source("R/lp_model.R")

# The state carried from one invasion to the next.
lp_state <- function(params, n, n_pp) list(params = params, n = n, n_pp = n_pp)

# Rebuild the params object for a given set of life histories.  The size grid
# is fixed (min_w = egg mass, max_w = 40 kg) and so does not move as species
# come and go, which lets abundance arrays be carried across rebuilds.
lp_rebuild <- function(w_max, mu_egg, z) {
    lp_params(w_max = w_max, mu_egg = mu_egg, z = z,
              species = seq_along(w_max))
}

# Draw one invader (Appendix B step 1).
lp_draw_invader <- function(assembly = LP_ASSEMBLY, randomise_A = FALSE) {
    list(w_max  = exp(runif(1, log(assembly$w_max_range[1]),
                            log(assembly$w_max_range[2]))),
         mu_egg = runif(1, assembly$mu_egg_range[1], assembly$mu_egg_range[2]),
         # Figs 5, 6: search rate and intrinsic mortality share a factor
         # z_i ~ N(1, 0.1), so faster feeders also die faster.
         z      = if (randomise_A) rnorm(1, 1, 0.1) else 1)
}

lp_assemble <- function(assembly = LP_ASSEMBLY, randomise_A = FALSE,
                        dt = 0.01, verbose = TRUE) {
    w_max <- numeric(0); mu_egg <- numeric(0); z <- numeric(0)
    n <- NULL; n_pp <- NULL
    attempts <- 0

    while (length(w_max) < assembly$n_species &&
           attempts < assembly$max_attempts) {
        attempts <- attempts + 1
        inv <- lp_draw_invader(assembly, randomise_A)

        w_max_new  <- c(w_max, inv$w_max)
        mu_egg_new <- c(mu_egg, inv$mu_egg)
        z_new      <- c(z, inv$z)
        p <- lp_rebuild(w_max_new, mu_egg_new, z_new)

        # Carry the resident spectra over; start the invader low on a power law.
        n_new <- initialN(p)
        if (!is.null(n)) n_new[seq_len(nrow(n)), ] <- unname(n)
        n_new[nrow(n_new), ] <- lp_invader_n(p, nrow(n_new), assembly)
        initialN(p) <- n_new
        if (!is.null(n_pp)) initialNResource(p) <- n_pp

        # Relax (Appendix B step 2)
        sim <- project(p, t_max = assembly$relax_years, dt = dt,
                       t_save = assembly$relax_years, progress_bar = FALSE)
        n_end    <- N(sim)[dim(N(sim))[1], , , drop = FALSE][1, , ]
        n_pp_end <- NResource(sim)[dim(NResource(sim))[1], ]
        if (is.null(dim(n_end))) n_end <- matrix(n_end, nrow = 1)

        # Extinction cull (Appendix B step 3)
        B <- drop(n_end %*% (w(p) * dw(p)))
        keep <- B >= assembly$extinct_B
        if (!any(keep)) {
            n <- NULL; n_pp <- n_pp_end
            w_max <- numeric(0); mu_egg <- numeric(0); z <- numeric(0)
        } else {
            w_max  <- w_max_new[keep]
            mu_egg <- mu_egg_new[keep]
            z      <- z_new[keep]
            n      <- n_end[keep, , drop = FALSE]
            n_pp   <- n_pp_end
        }
        if (verbose)
            message(sprintf("attempt %2d: invader w_max = %7.0f g -> %2d species",
                            attempts, inv$w_max, length(w_max)))
    }

    if (length(w_max) == 0)
        stop("assembly ended with no surviving species")

    # Order species by maximum body mass, as the paper numbers them (Fig. 2c).
    o <- order(w_max)
    p <- lp_rebuild(w_max[o], mu_egg[o], z[o])
    initialN(p) <- unname(n[o, , drop = FALSE])
    initialNResource(p) <- n_pp
    list(params = p, attempts = attempts)
}

# Relax an assembled ecosystem further, at the paper's time step, and return
# the params object holding the resulting quasi-equilibrium as its initial
# state.
lp_relax <- function(params, years = 50, dt = LP_NUMERICS$dt) {
    sim <- project(params, t_max = years, dt = dt, t_save = years,
                   progress_bar = FALSE)
    finalParams(sim)
}
