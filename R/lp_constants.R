# ---------------------------------------------------------------------------
# Law & Plank (2023) "Fishing for biodiversity by balanced harvesting",
# Fish and Fisheries 24, doi:10.1111/faf.12705.
#
# All parameter values from Tables 1 and 2 and Appendices B, C of the paper
# (preprint: doi:10.1101/2021.06.27.450047v3).  Symbols are named after the
# paper, with the paper's equation number given for each.
#
# Units throughout: grams, years, square metres of sea surface.
# ---------------------------------------------------------------------------

## Table 1: plankton ("resource" in mizer) -----------------------------------
LP_PLANKTON <- list(
    w_min   = 1e-10,  # g,  smallest cell           (w0 e^{x_0,0})
    w_max   = 1,      # g,  largest cell            (w0 e^{x_inf,0})
    r0      = 10,     # /yr at 1 g, intrinsic rate of increase      Eq. (A.12)
    rho     = 0.15,   # scaling of r with cell size                 Eq. (A.12)
    a0      = 2000,   # /m^2 at 1 mg, carrying capacity             Eq. (A.13)
    w_a     = 1e-3,   # g,  size at which a0 is quoted
    I0      = 2000,   # /m^2/yr at 1 mg, immigration rate           Eq. (A.14)
    w_I     = 1e-3,   # g,  size at which I0 is quoted
    lambda  = 2       # slope exponent of the isolated spectrum
)

## Table 2: fish -------------------------------------------------------------
LP_FISH <- list(
    w_egg     = 1e-3,   # g,  egg mass, same for all species
    beta      = 6.908,  # centre of the box feeding kernel, log(PPMR)  Eq. (A.4)
    sigma     = 1.535,  # width measure; kernel spans beta +/- 3 sigma
    A         = 37.5,   # m^2/yr/g^alpha, search-rate coefficient      Eq. (A.2)
    alpha_q   = 0.85,   # search-rate scaling exponent (mizer's q)
    K         = 0.1,    # conversion efficiency (mizer's alpha)
    w_L       = 0.1,    # g,  size at which larval mortality dies away Eq. (A.6)
    rho_L     = 5,      # sharpness of that transition                 Eq. (A.6)
    mu_b0     = 0.1,    # /yr, background mortality at egg size        Eq. (A.7)
    xi        = 0.15,   # size scaling of background mortality         Eq. (A.7)
    rho_m     = 15,     # sharpness of the maturation ogive            Eq. (A.8)
    rho_inf   = 0.2,    # approach to w_max                            Eq. (A.8)
    eps_R     = 0.2,    # reproductive efficiency                     Eq. (A.9)
    mu_b_form = "ratio",# reading of Eq. (A.7); see docs/mu_b.md
    w_mat_rat = 0.1,    # w_mat = w_max / 10                       (Appendix C)
    theta_ii  = 0.5,    # cannibalism weight                       (Appendix C)
    theta_ij  = 0.2     # between-species predation weight         (Appendix C)
)

## Appendix B: assembly and randomisation ------------------------------------
LP_ASSEMBLY <- list(
    w_max_range   = c(100, 40000),  # g,  log-uniform maximum body mass
    mu_egg_range  = c(28, 32),      # /yr, uniform larval death rate at egg size
    n_egg_init    = 0.002,          # /m^2, log-density of eggs of an invader
    invader_slope = 2,              # invader starts on u(x) ~ w^-(lambda-1)
    relax_years   = 50,             # yr of relaxation after each invasion
    extinct_B     = 2e-6,           # g/m^2, biomass below which a species is cut
    n_species     = 15,             # target richness
    max_attempts  = 40              # cap on invasion attempts
)

## Appendix C: numerics ------------------------------------------------------
LP_NUMERICS <- list(
    dx = 0.1,     # log body-mass step
    dt = 0.002    # yr, time step (Euler)
)

## Section 2.3: fishing ------------------------------------------------------
LP_FISHING <- list(
    w_f     = 400,   # g,  knife-edge entry size of the mixed-species fishery
    t_max   = 50,    # yr of fishing
    # Figs 3, 4: no randomisation of fishing intensity
    F_base  = 0.1,   # /yr                                          Eq. (2.7)
    c_P     = 1,     # m^2/g                                        Eq. (2.8)
    c_PB    = 0.25,  # dimensionless exploitation ratio             Eq. (2.9)
    # Fig. 6: doubled baseline, random intensity factor z' ~ U(0.5, 1.5)
    F_base6 = 0.2,
    c_P6    = 2,
    c_PB6   = 0.5,
    z_range = c(0.5, 1.5)
)
