# Why production is proportional to biomass

[Figure 2b of the overview](index.md) shows biomass and production rate
tightly correlated across species, over six decades. That relationship is the
premise the whole paper rests on: it is why $$F_i \propto P_i$$ makes fishing mortality
track abundance, and why $$F_i \propto P_i/B_i$$ barely differs from a constant
$$F$$. This page asks where it comes from, and then tries twice to break it.

Everything below is measured, not argued.

Four results:

1. $$P/B$$ is a community property, not a species property. It equals the
   biomass-weighted mean death rate, and mortality in this model is a function
   of body size that no species can escape. Abundance cancels exactly, so a
   species six decades rarer has the same turnover rate.
2. Making **egg mass** an independent life-history axis does not spread
   $$P/B$$. Assembly removes the variation first: egg mass is a fitness axis,
   and every surviving species is driven to the largest egg on offer.
3. Making **activity** an independent axis does work, but only once it is given
   a cost. Coupling a species' vulnerability to predation to its activity makes
   the axis exactly neutral, and then a genuine fast-slow continuum assembles —
   a ten-fold range in age at maturity, $$P/B$$ spread 3.5x wider in logs than
   the paper's. And $$B \sim P$$ *still* holds, with $$\alpha$$ at 0.973-1.003.

4. Run the harvesting rules on the result and BH<sub>P/B</sub> is no longer a
   near-duplicate of a fixed $$F$$ — it spans an order of magnitude and tracks
   turnover at $$r = 0.94$$ — but it still fails to protect rare species,
   because turnover is orthogonal to rarity. That is a stronger vindication of
   the paper's conclusion than its own ecosystems provide. BH<sub>B</sub>
   likewise still matches or beats BH<sub>P</sub>, on an assemblage built to
   give production something to do.

Reproduce with `Rscript run_egg_assembly.R <seed> <varied|fixed>`,
`Rscript run_activity_assembly.R <seed> <neutral|gradient>`,
`Rscript run_activity_harvest.R <seed>` and `Rscript run_activity_bhb.R <seed>`.

---

## The names and symbols used here

The model is Law & Plank's, written out in full in
[the overview](index.md). The one feature of it that this whole page turns on
is that **species differ in only two drawn quantities** — maximum body mass
$$w_{\max,i}$$ and larval death rate $$\mu_{\mathrm{egg},i}$$ — and are
identical in everything else, including egg mass. The experiments below add a
third axis and see what happens.

| | |
|---|---|
| `eco1` | the 15-species assemblage behind the paper's Figs 2-4, grown from seed 101 with the same search-rate coefficient $$A$$ for every species |
| `eco_r1`, `eco_r2`, `eco_r3` | three further assemblages (seeds 201-203) with the search rate randomised per species, $$A_i = z_i A$$, $$z_i \sim N(1, 0.1)$$ |
| `varied_301`, `fixed_301`, … | egg-mass runs, named *arm* underscore *seed*: `varied` draws an egg mass per species, `fixed` keeps the paper's common $$w_0$$ |
| `neutral_301`, `gradient_301`, … | activity runs: `neutral` gives the activity axis a predation cost, `gradient` is the paper's own treatment |
| $$B_i$$, $$P_i$$ | biomass and somatic production rate of species *i*, both measured over the harvested range $$w > w_f = 400$$ g (Eqs 2.3, 2.4) |
| $$\alpha$$ | the regression slope of $$\log B_i$$ on $$\log P_i$$ across the species of one assemblage. $$\alpha \approx 1$$ is the paper's premise, and every experiment below is a way of trying to move it |
| $$z_i$$ | activity factor: multiplies species *i*'s search rate **and** its intrinsic mortality, so faster feeders also die faster |
| $$\theta$$ | the interaction matrix — 0.5 on the diagonal (cannibalism), 0.2 off it. $$\theta_{ij}$$ is how readily *i* eats *j*, so column *j* is how readily *j* is eaten |
| $$R_0$$ | lifetime reproductive output of a rare invader; $$R_0 = 1$$ is the invasion threshold |
| `LP_FISH`, `vuln_exp` | the fish parameter list in [`R/lp_constants.R`](https://github.com/gustavdelius/balanced-harvest/blob/main/R/lp_constants.R), and the switch in it that couples vulnerability to activity. `vuln_exp = 0` is the paper |

## 1. Where the relationship comes from

### P/B is the mean mortality rate

For an unfished population in balance every gram produced is eventually lost,
so production equals losses:

$$
P_i = \int \mu_i(w)\, b_i(w)\, \mathrm{d}w
\qquad\Longrightarrow\qquad
\frac{P_i}{B_i} = \langle \mu_i \rangle_B ,
$$

the biomass-weighted mean death rate. On `eco1` the identity holds to 2-3%,
the residual being the quasi-equilibrium drift the paper acknowledges (the
worst offender is species 10, which is nearly extinct):

| species | 1 | 5 | 9 | 12 | 15 |
|---|---|---|---|---|---|
| $$P/B$$ | 0.996 | 0.704 | 0.682 | 0.611 | 0.568 |
| $$\langle\mu\rangle_B$$ | 1.019 | 0.728 | 0.702 | 0.628 | 0.581 |

So Figure 2b is really asking why every species has nearly the same average
mortality rate.

### Mortality and growth are functions of size, not of species

A fish's per-capita rates are set by where it sits in the community spectrum,
not by what it is. Predation mortality comes from the whole assemblage through
an almost uniform interaction matrix. Larval mortality acts only below
$$w_L = 0.1$$ g. And the background term, Eq. (A.7), depends on the species
only through how well it is feeding:

$$
\mu^{b}_i(w,t) = \mu_b^{(0)} \left(\frac{w}{w_0}\right)^{-\xi}
  \frac{\tilde g_i(w_0,t)}{\tilde g_i(w,t)},
  \qquad
  \mu_b^{(0)} = 0.1\ \mathrm{yr^{-1}},\quad \xi = 0.15,
$$

with $$\tilde g_i = K E_i/w$$ the mass-specific rate of food intake and
$$w_0 = 10^{-3}$$ g the egg mass. So it is a function of body size and of local
feeding conditions, and those are community properties: a fish of a given size
eats what the spectrum around it offers, whatever species it belongs to. (This
is the "ratio" reading of an equation the paper prints the other way up; see
[mu_b.md](mu_b.md).)

Measured across the 15 species of `eco1` at 0.2-30 g — past the larval term,
below every species' $$w_\mathrm{mat}$$:

| | max/min across species |
|---|---|
| mass-specific growth $$g(w)/w$$ | **1.14** |
| total mortality $$\mu(w)$$ | **1.27** |

against a spread in biomass of $$6\times10^5$$.

### Abundance cancels

Larval mortality sets the *amplitude* of $$b_i(w)$$, and amplitude divides out
of $$P/B$$ exactly. Nothing brings it back: there is no stock-recruitment
relationship, so rarity does not relax a recruitment brake, and there is no
satiation, so rarity does not buy better feeding. A rare species is too small a
part of the community to change the fields of growth and mortality it lives in.
Regressing $$\log(P/B)$$ on $$\log w_\mathrm{max}$$ and $$\log B$$, the
abundance term is 0.008 ± 0.004 (p = 0.10) — nothing, across six decades of it.

### What is left is w_max, and it is weak

$$\mu(w)$$ declines allometrically (log-log slope $$\approx -0.19$$ over
0.2 g - 20 kg), so a larger species puts more of its biomass where turnover is
slower:

$$
P_i/B_i \;\approx\; c\, z_i\, w_{\max,i}^{-0.14}
$$

| | eco1 | eco_r1 | eco_r2 | eco_r3 |
|---|---|---|---|---|
| $$P/B$$ spread | 1.76x | 1.67x | 2.00x | 1.68x |
| biomass spread | 5.8 dec | 3.0 | 3.8 | 7.8 |
| exponent on $$w_{\max}$$ | −0.142 | −0.140 | −0.161 | −0.069 |

Two decades of $$w_{\max}$$ buy at most a two-fold spread in $$P/B$$. On a
log-log plot that is a quarter of a decade of scatter against three to eight
decades of range — invisible, hence a line of slope 1. The residual is also
exactly what puts $$\alpha$$ just below 1 rather than at it: $$P/B$$ falls with
$$w_{\max}$$ and the large species are the rare ones, so $$P/B$$ is weakly
correlated with $$B$$ (0.46 in logs).

### Production is not where the biomass is

One decomposition matters for everything that follows. Splitting each species'
biomass and its losses by size band on `eco1`:

| | share of biomass | share of production |
|---|---|---|
| sub-gram fish | 0.4 - 2.3% | **18 - 59%** |
| above 400 g | 0 - 89% | 0 - 74% |

Sub-gram fish hold almost none of the biomass and generate up to 59% of the
production, the share falling monotonically as $$w_{\max}$$ rises. That
crossover *is* the $$w_{\max}^{-0.14}$$ relation, and it says the strongest
levers on $$P/B$$ are at the egg-and-larva end.

---

## 2. Egg mass: a lever that assembly removes

The paper gives every species the same egg mass $$w_0 = 10^{-3}$$ g (its
Table 2). Real egg sizes span decades and are the classic life-history axis
decoupled from adult size, so this is the obvious thing to vary.

The reason it ought to move $$P/B$$ is that larval mortality, Eq. (A.6), is a
reverse sigmoid in *absolute* body mass,

$$
\mu^{L}_i(w) = \frac{\mu_{\mathrm{egg},i}}{1 + (w/w_L)^{\rho_L}},
\qquad w_L = 0.1\ \mathrm{g},\quad \rho_L = 5,
$$

with $$\mu_{\mathrm{egg},i} \approx 30$$ yr⁻¹ — two orders of magnitude above
anything an adult experiences. A species with a 0.1 g egg is therefore born at
the *top* of that gauntlet, skipping the size range that generates a third to a
half of a small-egg species' production.

Egg mass was drawn log-uniform over $$[10^{-3}, 10^{-1}]$$ g, upwards from the
paper's value to $$w_L$$. Upwards because mizer takes the size grid's floor
from `min(species_params$w_min)`, so going below $$w_0$$ would move the grid —
and move it again whenever a species was added or culled, breaking the fixed
grid that assembly relies on. Egg masses are assigned onto the fixed grid after
construction and snapped to grid points.

Three seeds were assembled with eggs varied and three with the paper's common
$$w_0$$, so the comparison is not confounded by drawing a different assemblage.

### The axis collapses

Across all 27 surviving species in the varied runs:

- every survivor came from the **top 30%** of the drawn range
- median survivor sits at the **91st percentile**
- realised spread **0.61 of the 2 decades drawn**, piled against the cap
- KS against the uniform draw: D = 0.745, p = 1.9e-13

That is total directional selection toward the largest egg on offer, pressed
against the boundary — so the only thing stopping it going further is the range
limit, not anything in the model.

### P/B narrowed rather than widened

| | $$P/B$$ fold-range | sd($$\log_{10} P/B$$) | $$\alpha$$ | species |
|---|---|---|---|---|
| fixed eggs | 1.53, 1.95, 2.24 → **1.91** | **0.092** | 1.020 | 15, 15, 15 |
| varied eggs | 1.37, 1.62, 1.54 → **1.51** | **0.055** | 0.958 | 7, 9, 11 |

The strongest available lever, given two decades to work with, delivers *less*
$$P/B$$ variation than the paper's own parameterisation.

### Why: egg mass buys fitness, not turnover

Varying one species' egg mass in an otherwise fixed community, and re-relaxing:

| $$w_\mathrm{egg}$$ (g) | 0.001 | 0.005 | 0.020 | 0.049 | 0.090 |
|---|---|---|---|---|---|
| $$B_i$$ | 0.99 | 1.90 | 3.57 | 5.31 | **6.45** |
| $$P/B_i$$ | 0.671 | 0.562 | 0.620 | 0.760 | 0.925 |
| $$P/B$$, all others | 0.668 | 0.667 | 0.661 | 0.645 | 0.596 |

A 90-fold egg increase buys 6.5x the biomass and only 1.4x the $$P/B$$, and the
controlled slope $$\mathrm{d}\log(P/B)/\mathrm{d}\log w_\mathrm{egg} = 0.069$$
is four times weaker than the +0.27 measured across surviving species — most of
that cross-sectional association was confounding. The relationship is also
non-monotonic: $$P/B$$ falls at first, as the truncation argument predicts, then
reverses once the species becomes abundant enough to restructure the assemblage
around it (note the last row).

Egg mass is a *fitness* axis wearing a turnover axis's clothes. Assembly
climbed it to the ceiling, and the $$P/B$$ variation never materialised.

---

## 3. Activity: a lever that works, once it is given a cost

The model already has a fast-slow axis, and the paper already uses it: the
activity factor $$z_i$$ is what makes `eco_r1`-`eco_r3` differ from `eco1`.
It enters by scaling a species' search rate and its intrinsic mortality
together,

$$
A_i \to z_i A, \qquad
\mu_{\mathrm{egg},i} \to z_i \mu_{\mathrm{egg},i}, \qquad
\mu_b^{(0)} \to z_i \mu_b^{(0)},
$$

so that "faster feeders also die faster" — but it does **not** scale the
predation that other species impose on it, which is the column of $$\theta$$
belonging to that species.

That leak is fatal. Dropping a probe species into the assembled `eco1` at
negligible abundance and computing its lifetime reproductive output:

| $$z$$ | 0.35 | 0.60 | 1.00 | 1.70 | 2.86 |
|---|---|---|---|---|---|
| $$R_0$$, paper's treatment | 0.0020 | 0.061 | 0.408 | 1.326 | 2.623 |
| $$R_0$$, vulnerability coupled | 0.40824 | 0.40824 | 0.40824 | 0.40824 | 0.40824 |

$$
\frac{\mathrm{d}\log R_0}{\mathrm{d}\log z} = +3.33
\qquad\text{versus}\qquad
0 .
$$

Over an eight-fold range of $$z$$ the paper's treatment spans a factor of 1300
in fitness. The same steep ladder shows up in the published ecosystems: the
realised $$z$$ in `eco_r1..r3`, drawn from $$N(1, 0.1)$$, has mean 1.036
(p = 0.011).

The leak is small but sits in the wrong place. Biomass-weighted, predation is
only **9.5%** of total mortality (background $$\mu_b$$ 50.6%, larval 39.9%), so
the trade-off is already 90% complete. But at 1-10 g, where juvenile
survivorship is decided, predation is **93-95%** of mortality.

### The fix

Scale each species' vulnerability with its activity — multiply its *column* of
$$\theta$$ by $$z_i^{\,\phi}$$, so a species that forages harder is also more
exposed. This is the growth-predation-risk trade-off, not a fudge.

At $$\phi = 1$$ every rate a species experiences scales with $$z$$, so its life
is a pure time-rescaling. Survivorship $$\exp(-\int \mu/g\,\mathrm{d}w)$$ is
unchanged because both scale; reproductive output per unit time scales as
$$z$$; time spent per size interval as $$1/z$$. The two cancel identically, and
$$R_0$$ is invariant to five significant figures — there is no gradient left to
climb. It is set by `vuln_exp` in `LP_FISH`, default 0, which is the paper.

Because $$P/B = \langle\mu\rangle_B$$ and every mortality scales with $$z$$
while the biomass distribution keeps its shape, the prediction is
$$P/B \propto z$$ exactly.

### The neutral axis survives assembly intact

$$z$$ log-uniform over $$[0.35, 2.86]$$, matched pairs on three seeds:
`neutral` ($$\phi = 1$$) against `gradient` ($$\phi = 0$$, the paper's).

| | realised $$z$$ | median quantile of the drawn range | KS vs uniform |
|---|---|---|---|
| **neutral** | 0.35 - 2.75 | **0.474** | D = 0.116, **p = 0.55** |
| **gradient** | 1.49 - 2.85 | 0.903 | D = 0.706, **p = 4.7e-11** |

The neutral arm's survivors are statistically indistinguishable from the draw.
The gradient arm collapsed to the top, with the same signature egg mass gave.

### Richness recovers completely

| | species | attempts |
|---|---|---|
| neutral | **15, 15, 15** | 23, 24, 23 |
| gradient | 7, 6, 8 | 40, 40, 40 (cap) |

23-24 attempts is normal — the paper's own ecosystems took 23-33. So the
*fitness gradient* was what crippled richness, not the trait variation. That
also explains the egg-mass runs: those communities were not failing because
life-history variation destabilises coexistence, but because most invaders were
drawn onto a losing point of a steep ladder.

### P/B spreads, by the predicted amount

$$
\frac{\mathrm{d}\log(P/B)}{\mathrm{d}\log z} = 0.979 \pm 0.016
\qquad (\text{predicted } 1.0,\ R^2 = 0.99)
$$

| | $$P/B$$ fold-range | sd($$\log_{10} P/B$$) |
|---|---|---|
| the paper's four ecosystems | 1.67 - 2.00 | 0.058 - 0.078 |
| **neutral** | **5.23 - 7.75** | **0.237 - 0.255** |

And it is a real fast-slow continuum (`neutral_301`):

| | slowest (sp 8) | fastest (sp 2) |
|---|---|---|
| $$z$$ | 0.35 | 2.32 |
| age at maturity | **8.3 yr** | **0.85 yr** |
| biomass turnover time | 4.22 yr | 0.54 yr |

A ten-fold range in age at maturity, coexisting at 15 species. And
$$\mathrm{cor}(\log z, \log w_{\max}) = -0.085$$: the axis is independent of
body size, so the assemblage has fast small species and slow large ones rather
than a repackaging of $$w_{\max}$$.

### But B ~ P does not break

$$\alpha$$ comes back at **1.003, 0.994, 0.973** — if anything closer to 1 than
the paper's own 0.977-0.998.

The reason is arithmetic. $$\alpha$$ is the regression slope of $$\log B$$ on
$$\log P$$, so it is set by the ratio of the two variances, and biomass still
spans 4.9 decades while $$P/B$$ now spans 0.8. Even a six-fold $$P/B$$ spread
is small next to that. To move $$\alpha$$ appreciably you would need $$P/B$$
varying over a range comparable to biomass — orders of magnitude, not factors
of six. What did change is the scatter: $$\mathrm{cor}(\log B, \log P)$$ falls
from 0.997-1.000 to 0.970-0.991.

So $$B \sim P$$ survives a genuine fast-slow continuum. It is more robust than
the mechanism above suggests, and this strengthens the paper's premise.

---

## 4. The three rules on a fast-slow ecosystem

The middle row of [Figure 3](index.md) shows BH<sub>P/B</sub> allocating
nearly the same fishing mortality to every species, and that is the paper's
argument for why a constant exploitation ratio is barely distinguishable from
a constant $$F$$. The argument rests on $$P/B$$ being near-constant, which is a
property of the paper's parameterisation rather than of size-spectrum
ecosystems in general — so it should not survive here.

All three rules were calibrated to the same total yield at year 50 against
$$F = 0.1$$, exactly as the paper does it, on each of the three neutral
ecosystems (`Rscript run_activity_harvest.R <seed>`).

### BH_P/B stops being a near-duplicate of a fixed F

| | spread of $$F_i$$ across species | $$\mathrm{cor}(\log F_i, \log z_i)$$ |
|---|---|---|
| fixed | 1.0x by construction | — |
| BH<sub>P</sub> | 4.0e3 - 2.0e6 x | −0.18, −0.02, +0.13 |
| BH<sub>P/B</sub> | **5.7x, 5.8x, 9.1x** | **+0.96, +0.93, +0.94** |

BH<sub>P/B</sub> is now a genuinely distinct rule that fishes fast species six
to nine times harder than slow ones, at a correlation with activity of 0.94.
That is textbook single-stock advice: high-turnover stocks sustain higher
$$F$$. BH<sub>P</sub>, by contrast, is essentially uncorrelated with turnover —
it tracks abundance, at $$\mathrm{cor}(\log F_i, \log B_i) = 0.94$$ to
$$0.98$$.

### And it still fails, for a reason that makes the paper's case stronger

Speed and rarity are nearly orthogonal in these ecosystems:
$$\mathrm{cor}(\log z, \log B) = -0.32, -0.06, -0.12$$. So a rule that
allocates effort along the turnover axis is allocating along an axis that
carries almost no information about which species are at risk.

The worst-affected species under each rule says it plainly (rank 1 = rarest):

| | seed 301 | seed 302 | seed 303 |
|---|---|---|---|
| fixed | B rank 1 | 12 | 2 |
| BH<sub>P</sub> | **15** | **15** | **13** |
| BH<sub>P/B</sub> | 4 | 1 | 4 |

This is a **stronger** test of the paper's conclusion than the paper's own.
There, BH<sub>P/B</sub> could be dismissed as barely differing from a constant
$$F$$, so its poor showing proved little about the rule itself. Here it is
non-degenerate, spans an order of magnitude, and implements orthodox fisheries
advice — and it still leaves the rare species exposed, because turnover is not
what rarity is made of. Only tracking abundance protects biodiversity.

### The headline metric is misleading, and this is where it shows

Under BH<sub>P</sub> the worst-affected species is the *most abundant* one, by
design: effort is deliberately concentrated there. So "worst species relative
to control", the metric used in [the overview](index.md) and by the paper,
measures harm to a common species under one rule and harm to a rare one under
another, and the comparison between them means less than it appears to.

On that metric BH<sub>P</sub>'s advantage over fixed $$F$$ looks like
0.95x, 1.40x, 1.60x. Conditioned on rarity it looks quite different:

| advantage of BH<sub>P</sub> over fixed $$F$$ | seed 301 | seed 302 | seed 303 |
|---|---|---|---|
| worst species overall (the headline metric) | 0.95x | 1.40x | 1.60x |
| worst of the rarest third | **5.4x** | **10.8x** | **16.9x** |
| geometric mean over all 15 species | 1.38x | 1.83x | 1.64x |

Under BH<sub>P</sub> the rarest species ends up *above* its unfished
trajectory — 5.8 to 6.5x — because the rule barely fishes it while fishing its
competitors and predators hard.

Two caveats on that last number. In seeds 301 and 302 the single rarest species
sits at $$2.9\times10^{-5}$$ and $$5.6\times10^{-7}$$ g m⁻², at or below the
$$2\times10^{-6}$$ extinction threshold used during assembly, so the 18-23x
ratios there are measured on species that are already effectively gone. (The
threshold is applied during assembly only, and the final relaxation lets species
drift below it; `eco1` has the same property.) Only seed 303's rarest species,
at $$2.7\times10^{-4}$$ g m⁻², is comfortably above it — and there the
advantage is 17.7x on a species that genuinely persists. The rarest-third
column is the more robust statement.

### BH_P against BH_B: production still does no work

[robustness-results.md](robustness-results.md) finds that $$F \propto B_i$$
beats $$F \propto P_i$$ in 13 of 13 ecosystems. That could be read as an
artefact of $$P/B$$ being flat: if $$P_i = k B_i$$ with $$k$$ common to every
species, the two allocations are the same rule up to the calibration constant,
and any difference between them is noise. These ecosystems test that reading,
because here the two rules genuinely come apart:

$$
F_i(\mathrm{BH}_P) \propto P_i = (P_i/B_i)\, B_i \propto z_i B_i,
\qquad
F_i(\mathrm{BH}_B) \propto B_i ,
$$

so BH<sub>P</sub> is BH<sub>B</sub> tilted by the activity factor. Regressing
$$\log(F_i^{\mathrm{BH}_P} / F_i^{\mathrm{BH}_B})$$ on $$\log z_i$$ gives a
slope of **0.968 ± 0.080, 0.888 ± 0.105, 1.057 ± 0.110** against the 1.0 the
algebra demands, at $$R^2 = 0.86$$ to $$0.92$$. The tilt is real and it is
exactly $$z$$.

It does not help. A fourth rule, calibrated to the same yield target as the
other three:

| BH<sub>P</sub> relative to BH<sub>B</sub> (>1 favours BH<sub>P</sub>) | seed 301 | seed 302 | seed 303 |
|---|---|---|---|
| rarest species | 0.86x | **0.68x** | 0.74x |
| worst of the rarest third | 0.98x | 0.87x | 0.96x |
| geometric mean over 15 species | 1.00x | 0.99x | 1.01x |

BH<sub>B</sub> is the equal of BH<sub>P</sub> on the community as a whole and
better on the species that matter most, by 16% to 47% on the rarest one. So
multiplying the allocation by $$P/B$$ does not merely fail to help when $$P/B$$
is flat — it actively costs a little when $$P/B$$ varies by 5 to 8x, and the
metrics that weight rarity most are where the cost shows. The extra factor is
orthogonal to rarity, so it can only blur an allocation that should track
abundance alone.

The one column where BH<sub>P</sub> leads is "worst species overall", by 1.5x
to 1.9x — and that is the misleading metric again. Under both rules the
worst-affected species is the *most abundant* one (abundance rank 15, 15, 13),
and BH<sub>B</sub> presses harder on it precisely because it tracks abundance
without the tilt. That is the rule working as designed, not failing.

This strengthens the "production is doing no work" finding rather than
qualifying it: the result survives an assemblage built specifically to give
production something to do.

## Numerics

The fast end of the activity range runs about 2.9x the paper's rates, so the
time step was cut to match. With $$\phi = 1$$ a community in which every
$$z = Z$$ *is* the $$z = 1$$ community running $$Z$$ times faster, which makes
an exact convergence test available. Biomass after 20 years, relative to
$$\Delta t = 0.001$$:

| | $$\Delta t = 0.01$$ | $$\Delta t = 0.004$$ |
|---|---|---|
| $$Z = 1$$ (paper's speed) | 0.10% | 0.03% |
| $$Z = 2.86$$ (fast end) | 0.67% | 0.23% |

Assembly therefore uses $$\Delta t = 0.004$$ and the final relaxation
$$\Delta t = 0.001$$, which keeps the fast end within about twice the
discretisation error the paper's own $$\Delta t = 0.01$$ carries at $$Z = 1$$.

## Caveats

- Three seeds per arm is a small sample, and the gradient arms stopped at the
  40-attempt cap rather than at an equilibrium richness, so the richness
  comparison is a lower bound on the gap rather than a measured asymptote.
- Exact neutrality is degenerate: zero gradient means nothing selects *for*
  coexistence along the axis either, so species may drift along it. The
  realised distributions match the draw after assembly, but these runs do not
  establish that the range is stable over much longer relaxations. A $$\phi$$
  slightly below 1 would give weak stabilisation.
- Cannibalism breaks the rescaling slightly at realistic abundance: mortality
  on $$i$$ from predator $$j$$ scales as $$z_i$$ (column) times $$z_j$$ (search
  volume), so for $$j = i$$ it scales as $$z^2$$. Negligible for an invading
  probe, a small deviation once the species is abundant.
- The egg-mass range runs upward from the paper's value rather than straddling
  it, for the grid reason given above, so every species in those runs has at
  least the paper's egg mass.
- `varied_301` has a biomass span of 0.96 decades and so fails Appendix B's
  fourth selection criterion, that biomass span roughly four orders of
  magnitude across species. It is not a valid Law & Plank ecosystem and should
  not be used downstream.
