# Project notes: building the Law & Plank model

## A box feeding kernel is not normalised in mizer

`box_pred_kernel()` returns 1 between `ppmr_min` and `ppmr_max`. Law & Plank's
Eq. (A.4) has height `1/(6 sigma)` so that the kernel integrates to 1 over log
prey mass. Fold that into `gamma`:

```r
gamma = A / (6 * sigma)          # A = 37.5, sigma = 1.535
ppmr_min = exp(beta - 3*sigma)   # 10
ppmr_max = exp(beta + 3*sigma)   # 1e5
```

Check it by comparing `getEncounter()` against a hand-coded Eq. (A.2); see
`tests/test_model.R`.

## Reaching an arbitrary psi through species parameters

`psi = maturity * repro_prop` with
`maturity = (1 + (w/w_mat)^-U)^-1`, `U = log(3)/log(w_mat/w_mat25)`, and
`repro_prop = (w/w_repro_max)^(m - n)`. So Eq. (A.8)'s `rho_m = 15` and
`rho_inf = 0.2` are reached with

```r
w_mat25 = w_mat / 3^(1/15)       # gives U = 15
n = 2/3, m = 2/3 + 0.2           # gives (w/w_max)^0.2
```

No need to pass a `repro_prop` array. Note mizer rounds `maturity` below 1e-8
down to zero, so compare against the analytic form absolutely, not relatively.

## mizer's RDI carries a sex ratio the paper does not

`mizerRDI()` multiplies by 0.5 for a 50:50 sex ratio. Eq. (A.9) has no such
factor, so `erepro = 2 * eps_R`.

## Keep the size grid fixed during sequential assembly

Build every intermediate params object on the same grid — `min_w` = egg mass,
`max_w` = the largest possible `w_max`, `min_w_pp` = smallest cell — and let
species carry their own smaller `w_max`. Abundance arrays are then compatible
across rebuilds, so a species can be added or culled by `rbind`/subsetting
rather than by `addSpecies()`, which would re-derive defaults you have
deliberately overridden.
