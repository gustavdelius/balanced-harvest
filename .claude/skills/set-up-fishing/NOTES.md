# Project notes: fishing

## Balanced harvesting is a rate function, not a gear

This project implements balanced harvesting in the sense of Law & Plank (2023),
where fishing mortality is **independent of body size** above a single shared
entry mass w_f, and the balancing happens **across species**:

    F_i(t) = c_P P_i(t)          (BH_P)
    F_i(t) = c_P/B P_i(t)/B_i(t) (BH_P/B)

The selectivity side is therefore trivial — `knife_edge` with an explicit
`knife_edge_size`. Note that mizer's default `knife_edge_size` is `w_mat`,
which is *per species*; a single shared entry mass must be set explicitly or
you silently get a maturity-based fishery instead.

## mizerRates computes e_growth before f_mort

`mizerRates()` fills `r$e_growth` before it calls `FMort`, and passes it in as
the `e_growth` argument. So a custom FMort can compute a species' somatic
production rate from the state it is handed, with no extra work:

    P_i = sum_w e_growth[i, w] * n[i, w] * dw[w]   over the harvested range

This is what makes a state-dependent, adaptive fishing rule possible at all —
`project()`'s `effort` argument must be prescribed in advance and cannot see
the evolving state. See `bh_model.R:balancedFMort`.

## Don't set the selectivity array by hand for this

An earlier attempt built a productivity-proportional selectivity array and
pushed it in with `selectivity(params) <- ...`. That works but **freezes** the
array (mizer marks it `"set manually"` and stops rebuilding it from
`gear_params`), and it is static through a projection. For a feedback rule,
`setRateFunction(params, "FMort", ...)` is the right mechanism.

## Watch the units of c_P

c_P has dimensions area/mass, so its numerical value is tied to the model's
biomass units and is not portable. Law & Plank's c_P = 1 m^2 g^-1 corresponds
to ~1.3e-12 in NS_params, where biomass is in grams over the whole North Sea.
Always calibrate it inside the model at hand.

## w_f will not land on a grid point

`w_f = 400 g` is not on the log grid, so the fishery starts at the first grid
point at or above it. Use that same index for the selectivity, for `B_i` and
for `P_i` (`R/lp_harvest.R:lp_wf_idx`), otherwise `Y_i = F_i B_i`, Eq. (2.6),
stops holding exactly.

## getYield() ignores a custom FMort

`getYield()` recomputes fishing mortality per gear from effort x catchability x
selectivity, so it does not see a `FMort` registered with `setRateFunction()`.
With a state-dependent rule, compute the yield yourself as `F_i(t) B_i(t)`,
which is what Eq. (2.6) says anyway. See `R/lp_harvest.R:lp_yield`.
