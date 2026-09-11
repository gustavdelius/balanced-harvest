# Project notes: extending mizer

Findings from implementing Law & Plank (2023); see `docs/mapping.md`.

## The Mort rate function is not handed the encounter rate

`mizerRates()` computes `encounter` first, but only passes `f_mort` and
`pred_mort` on to `Mort`. A mortality term that depends on the current food
intake — here the background mortality of Eq. (A.7) — therefore has to
recompute it, one extra convolution per time step:

```r
lpMort <- function(params, n, n_pp, n_other, t = 0, f_mort, pred_mort, ...) {
    encounter <- mizerEncounter(params, n = n, n_pp = n_pp, n_other = n_other, t = t)
    pred_mort + params@mu_b + f_mort + lp_background_mort(params, encounter)
}
```

Caching the encounter between the two calls is tempting and unsafe: `getMort()`
can be called on its own, and the cache would then be stale. Precompute the
*static* factors instead (stored in `other_params`) and keep the convolution.

By contrast `FMort` *is* passed `e_growth`, which is what makes a fishing rule
that responds to the current production rate possible at all.

## Turning off satiation and metabolism

The paper has no maximum intake rate and no metabolic cost, so
`E = alpha * encounter`. Setting `h = Inf` gives a zero feeding level, but
registering a `FeedingLevel` function that returns `encounter * 0` avoids any
`Inf/Inf` risk. Metabolism goes away with `ks = 0, k = 0`.

## mizer's predation diffusion omits the psi factor

`mizerDiffusion()` returns `(1 - f) * search_vol * alpha^2 * I_d`, without the
fraction of food that goes to somatic growth. Law & Plank's Eq. (A.1f) carries
`epsilon_i = 1 - psi`, so multiply it back in. Remember
`use_predation_diffusion(params) <- TRUE` — it is off by default.

## Custom resource dynamics: prefer an exact step to Euler

The plankton logistic has `r` up to ~320/yr at 1e-10 g, so explicit Euler needs
`dt < 0.006`. `dn/dt = I + (r-d) n - (r/a) n^2` is a Riccati equation with
coefficients constant over a step, so it has a closed-form solution:

```r
k <- r/a; disc <- sqrt((r-d)^2 + 4*k*I)
n_hi <- ((r-d) + disc)/(2*k); n_lo <- ((r-d) - disc)/(2*k)
v <- (n_hi - n)/(n - n_lo) * exp(-disc*dt)
n_new <- (n_hi + n_lo*v)/(1 + v)
```

Written this way `v >= 0` decays, so there is no overflow and no stability
limit. Same fixed points as the ODE. (mizer's own `resource_logistic()` does
the same thing for the case without immigration.)

## The paper's log-densities are not mizer's densities

Law & Plank's state variable is `u(x) = w phi(w)`, a density per unit *log*
mass. Carrying capacity and immigration quoted for `u` must be divided by `w`
before they go into `resource_capacity` and the immigration vector. The check
that catches an error here is the total primary production, `int r(w) w phi dw`,
which has to come out near 4000 g/m2/yr.
