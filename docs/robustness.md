# Which arbitrary choices could change Law & Plank's conclusion, and how to find out

Status: the design document. Results as they come in are in
[robustness-results.md](robustness-results.md), which also corrects two claims
made here (section 3.1 on the calibration criterion, and the role of `I_0` in
section 3.6). It is written against
the reimplementation in this repository, and names the code hooks each
experiment would need.

## 1. The claim, stated precisely

> At matched total yield, setting fishing mortality in proportion to somatic
> production rate (`F_i = c_P P_i`) leaves rare species roughly where they were
> after 50 years, whereas a fixed `F` or a fixed exploitation ratio
> (`F_i = c_P/B P_i/B_i`) puts several of them on an exponential path to
> extinction.

The paper's stated mechanism is a feedback: because `B_i ~ P_i^alpha` with
`alpha` near 1, `F_i = c_P P_i` makes fishing mortality roughly proportional to
biomass, so it falls away as a stock falls. Fixed `F` has no such feedback, and
`P_i/B_i` changes little as a stock falls, so `BH_P/B` has almost none either.

## 2. The causal chain, and where it can be attacked

| # | Link | What it needs to be true | Attack |
|---|---|---|---|
| i | `B_i ~ P_i^alpha`, `alpha ~ 1` | enough adult mortality that biomass does not pile up at `w_max`; a sensible size range over which `P` is measured | Eq. (A.7); `w_f`; the P/B measurement range |
| ii | `F_i ∝ P_i` therefore `F_i ∝ B_i` | mass-specific production roughly density-independent | food-limitation feedback; `theta` |
| iii | that feedback is what saves rare species | rare species decline *because of fishing* | baseline drift; predation release |
| iv | the other two rules have no comparable feedback | no other density dependence in the model | absence of a stock-recruitment relationship |
| v | "at matched total yield" | the matching criterion is not itself doing the work | terminal-yield calibration; fishing intensity |

Links iii, iv and v are the ones the paper does not examine, and they are the
ones most likely to matter.

## 3. Choices that could change the conclusion

### 3.1 Matching total yield *after 50 years* (Section 3)

The three constants are calibrated so that all three regimes give the same total
yield **in year 50**. By year 50 the fixed-`F` and `BH_P/B` ecosystems have been
substantially depleted, so their terminal yield is low; matching it forces
`c_P` down and makes `BH_P` fish *gently throughout*. Matching year-0 yield, or
yield integrated over the 50 years, would each give a different `c_P` — plausibly
a much larger one.

This is the single most load-bearing arbitrary choice in the paper, because the
whole comparison is "at equal yield" and there are at least four defensible
meanings of that phrase. It is not discussed.

**Prediction.** `c_P` calibrated on cumulative yield will exceed `c_P` calibrated
on terminal yield, perhaps substantially, and the advantage of `BH_P` will
shrink. Whether it *reverses* is the question.

### 3.2 No compensatory recruitment (Appendix B)

The model deliberately has no stock-recruitment relationship: recruitment is
strictly proportional to adult reproductive output, `RDD = RDI`. The paper
presents this as a virtue — density dependence should emerge rather than be
imposed — and it is a defensible modelling stance.

But it is also precisely the assumption that makes exponential decline possible.
With linear recruitment, a species whose adult mortality exceeds what its
lifetime reproductive output can sustain declines geometrically, forever. Add
even moderate compensation (a Beverton–Holt `R_max`, i.e. a non-zero
`reproduction_level`) and every species acquires a density-dependent brake at
low biomass — under *all three* rules. The `BH_P` feedback would then be
supplying protection that the population already has.

Most published multispecies size-spectrum models (including mizer's North Sea
model) do impose such a relationship. So this is not an exotic alternative; it
is the mainstream one.

The authors are aware that something has to take the place of an imposed
stock-recruitment relationship, and nominate a specific mechanism for it — larval
competition for plankton (Appendix B, device 1). Section 3.6 shows, by
measurement, that at the published parameter values that mechanism is barely
operating.

**Prediction.** As `reproduction_level` rises from 0 towards 0.5, the fixed-`F`
collapses turn into new lower equilibria, and the gap between the three rules
narrows sharply. If it closes, the paper's result is a statement about
recruitment dynamics as much as about balanced harvesting.

### 3.3 Which species are declining anyway (Appendix B, Section 3)

The unexploited ecosystem is explicitly *not* at equilibrium: "some species'
biomasses were still changing slowly after 50 years, so the state at this time is
best thought of as a quasi-equilibrium". Sequential assembly with an extinction
threshold leaves some species close to exclusion, and some of the "rare" species
may have been on their way out before any fishing started.

Every figure compares year 50 *fished* against year 0 *unfished*, so baseline
drift and fishing effect are confounded. The fix is trivial — run an unfished
control over the same 50 years and compare fished against control rather than
against year 0 — and this repository already does it, but the paper does not.

**Prediction.** Some of the apparent fixed-`F` collapse is drift. The ranking of
the three rules probably survives, but the magnitudes will shrink, and possibly
by a lot for the rarest species.

### 3.4 Is it production, or just density dependence?

`BH_P` gives `Y_i = c_P P_i B_i`, and since `P_i ~ B_i`, `Y_i ~ B_i^2`. A rule
`F_i = c_B B_i` gives exactly the same thing without ever mentioning production.
If the two behave identically, then the paper's result is about making fishing
mortality density-dependent, and "production" is a convenient proxy for biomass
rather than the active ingredient — which matters, because biomass is far easier
to estimate than production.

More sharply, the three rules are nearly the same one-parameter family:

```
F_i = c X_i^theta      with  theta = 0 (fixed F),  theta ~ 1 (BH_P),  theta <~ 0 (BH_P/B)
```

`BH_P/B` sits at `theta` slightly *negative*, because as a stock thins, competition
for food eases, `P_i/B_i` rises, and fishing mortality goes *up* — the paper says
as much ("The ratio might even increase if low species biomass reduced
competition for food"). Sweeping `theta` continuously turns a three-way
categorical comparison into a quantitative question: how much density dependence
does a harvest control rule need?

### 3.5 Adaptive feedback, or just a better allocation?

The paper's policy claim is specifically about adaptation over time: fishing
"operating adaptively to follow species' production rates over time, contains a
feedback that would help to protect species from overfishing in the presence of
uncertainty".

But `BH_P` does two things at once. It *allocates* effort across species at the
outset (rare species get low `F` from year 0), and it *updates* that allocation
as the ecosystem changes. These can be separated, and the paper does not separate
them. The decisive experiment is a frozen control: set `F_i = c_P P_i(0)`,
constant in time. If frozen `BH_P` does as well as adaptive `BH_P`, the benefit is
the allocation, and the feedback story is decoration.

There is also a second, indirect route the paper mentions only in passing:
"we found cases where rare species gained more through release from predation
under `BH_P` than they lost from fishing (results not shown)". If the protection
of rare species is mostly trophic release — because `BH_P` fishes the abundant
large predators hardest — then it depends on the model's predation structure
rather than on the harvest rule's feedback, and would not transfer to an
assemblage with a different trophic arrangement. This is decomposable: compare
each rare species' direct fishing loss against its change in predation mortality.

### 3.6 The plankton is almost impossible to eat down, so larval competition barely operates

Appendix B nominates larval competition as the first of the three devices that
make coexistence possible without an imposed stock-recruitment relationship:

> When larval density is low relative to plankton food, the abundance of food
> allows fast body growth through the vulnerable larval stage, making the risk of
> death relatively small. Conversely, when larval density is relatively high,
> growth is slower and the accumulated risk of death in the larval stage becomes
> greater. This has a strong stabilizing effect [...]

That feedback exists only to the extent that larvae actually deplete their food.
They do not. Measured on the assembled 15-species ecosystem `eco1`:

| plankton size | `r(w)` /yr | grazing mortality /yr | resource level `n/a` |
|---|---|---|---|
| 1e-8 g | 158 | 0.04 | **1.006** |
| 1e-6 g | 79 | 0.69 | **1.004** |
| 1e-5 g | 56 | 0.73 | **1.005** |
| 1e-4 g | 40 | 0.87 | **1.003** |
| 1e-3 g | 28 | 1.17 | **0.994** |
| 0.01 g | 20 | 9.73 | 0.596 |
| 0.1 g | 14 | 10.09 | 0.445 |
| 1 g | 10 | 9.91 | 0.321 |

A 1 mg larva eats prey between 1e-8 and 1e-4 g. Over that entire range the
plankton sits **at its carrying capacity**, despite fifteen fish species grazing
it. Only above about 0.01 g, where fish of 0.1 g and up feed, does grazing bite —
and there it bites hard, pulling the top of the plankton spectrum down to a third
of capacity.

The reason is structural, and follows from Eq. (A.11) in two lines. At steady
state, dividing through by `a(x)` and writing `L = n/a` and `iota = I/a`:

```
iota + r L (1 - L) - d L = 0        =>    r L^2 - (r - d) L - iota = 0
```

With `r = 10 w^-0.15`, that is 40-220/yr across the larval prey range, against
a measured grazing mortality of 0.4-0.9/yr. Solving for `L` at `r = 79/yr`,
`d = 0.7/yr`:

| `iota` (per year) | 1 (published) | 0.1 | 0.01 | 0 |
|---|---|---|---|---|
| resource level `L` | 1.004 | 0.992 | 0.991 | 0.991 |

**The immigration term is not what protects the plankton.** Switching it off
entirely moves the resource level from 1.004 to 0.991. What keeps `L` pinned at
1 is that the plankton regenerate two orders of magnitude faster than they are
grazed. (An earlier draft of this section named `I_0 = a_0` as the decisive
choice; that was wrong. It becomes decisive only in combination with a slow
plankton — at `r = 0.79/yr` the same table reads 1.183, 0.417, 0.183, 0.114.)

So the parameter to reach for is `r_0`, and `iota` matters only once `r_0` is
already low. Both are set without justification: `r_0 = 10` is anchored to a
cell-division scaling exponent, not to a rate, and `I_0` is simply set equal to
`a_0` with no comment beyond the statement that `I(x)` "was assumed to scale
with body size in the same way as the carrying capacity".

The effect on growth is direct. Multiplying the whole fish spectrum by 2 and
letting the plankton re-equilibrate changes growth rate by:

| fish body mass | d log(growth) / d log(fish abundance) |
|---|---|
| 0.001 g (egg) | **-0.013** |
| 0.01 g | -0.021 |
| 0.1 g | -0.093 |
| 1 g | -0.271 |
| 10 g | -0.524 |
| 100 g | **-0.908** |

Doubling the fish nearly halves the growth rate of a 100 g fish and slows an egg
by 1.3%. Competition for food is a real and very strong force in this model —
but only from about 0.1 g upwards, two orders of magnitude above the egg.

Three consequences.

1. **Egg-to-recruit survival is nearly independent of larval density *through
   growth*.** Since the model also has no satiation (`f == 0`), larval growth rate
   is essentially a fixed constant, so time spent in the vulnerable stage is
   fixed, so the growth-mediated part of the larval feedback contributes almost
   nothing. What density dependence remains at the larval stage comes through
   *predation* on eggs and larvae by larger fish, which is a different mechanism
   from the one Appendix B describes.

2. **It strengthens the worry in §3.2 rather than relieving it.** The model lacks
   both an imposed stock-recruitment relationship and (at these parameters) the
   emergent larval mechanism offered in its place. The recruitment brake is
   therefore weaker than the paper's own account of the model implies, which makes
   exponential decline under fixed `F` easier and makes the `BH_P` feedback look
   more necessary than it would in a model with a working larval feedback.

3. **The absence of satiation cuts the other way, and loses.** With no maximum
   intake rate, growth is strictly proportional to food, so a 10% drop in food is a
   10% drop in growth; in a standard mizer model at a feeding level of 0.6 it would
   be about 4%. Law & Plank's choice therefore *amplifies* whatever competition
   exists — and the elasticity at egg size is still only -0.009.

**An aside from the same measurement.** Averaged over the fifteen species,
plankton supplies 100% of intake at egg size and still 98% at 1 kg. By mass
flow, every fish in this model is a planktivore. That is a consequence of the
box kernel spanning four decades of prey size — the abundant small end dominates
the integral — together with the plankton extending up to 1 g. Species are
strongly coupled in this model, but through *predation mortality* rather than
through competition for fish prey, and a 1 kg fish sits at a trophic level near
2 rather than at the top of a chain. Anything that depends on the length of the
food chain (trophic cascades, indirect effects of fishing one species on
another) should be expected to be sensitive to the kernel width and to the
plankton cutoff of §4.4.

**What to vary.** `r_0` and `rho` (how fast plankton regenerate), `I_0/a_0` (the
floor that grazing cannot push through), `a_0` (standing stock), and `lambda`
(which shifts food between small and large plankton). Note that `r` and `I/a`
control the resource *level* directly, whereas `a_0` acts indirectly — raising it
raises fish biomass and hence grazing pressure too — so `a_0` is closer to a
nuisance scale parameter than the other two. Lowering `I_0/a_0` by one or two
orders of magnitude, or lowering `r_0`, should switch the larval feedback on and
is the cleanest test of whether the paper's stated coexistence mechanism, once
actually operating, changes the comparison between rules.

## 4. Choices that could change the magnitude

### 4.1 `w_f = 400 g`, a single knife edge for every species

The paper flags this as a simplification ("we picked a single illustrative
pattern of fishing [...] there are clearly many other possibilities"), but does
not test it, and it is more consequential than it looks. With `w_mat = w_max/10`
for every species:

* a 1 kg species matures at 100 g and is fished only as an adult;
* a 40 kg species matures at 4 kg and is fished for a decade before it breeds.

So a fixed `w_f` combined with a fixed `w_mat/w_max` ratio *systematically*
subjects large species to recruitment overfishing and shields small ones. Some of
"large species are vulnerable, small ones are not" is built in by this choice
rather than discovered. Species with `w_max < 400 g` are not fished at all.

### 4.2 Measuring `P` and `B` over the harvested range only

Justified on data-availability grounds ("reliable information is most likely to
be available over this range"), but it makes `P_i` a strange quantity: it excludes
juvenile production, which is most of a species' somatic production, and for a
species whose `w_max` is only a little above `w_f` it is dominated by the
boundary influx term at `w_f`. That makes `P_i` — and hence `F_i` under `BH_P` —
very sensitive to where `w_f` falls in a species' life cycle.

### 4.3 The interaction matrix (Appendix C)

`theta_ii = 0.5`, `theta_ij = 0.2`, flat, chosen to give diagonal dominance
because that "is known to promote coexistence". Two issues. First, the 2.5:1
ratio is not measured, and coexistence — hence which species are rare — depends
on it. Second, strong cannibalism is itself a density-dependent brake, so it may
be supplying some of the stability attributed to `BH_P`. Real matrices (Blanchard
et al. 2014; Spence et al. 2021) are heterogeneous, and with a flat off-diagonal
all 15 species are ecologically interchangeable apart from `w_max`.

### 4.4 Plankton capped at 1 g

(See also §3.6, which is about the same spectrum but a different property of it.)

Admitted to be artificial, "to compensate for the absence of multicellular
zooplankton". With a feeding kernel spanning prey from 1/100000 to 1/10 of body
mass, every fish above about 10 g must get part of its diet from fish. Raising
the cap weakens fish–fish coupling and therefore weakens both the trophic
cascades and the predation-release effect of §3.5.

### 4.5 Fishing intensity

The paper already notes that adaptive fishing "worked better at lower, than at
higher fishing intensities". The baseline is doubled from 0.1 to 0.2/yr between
Fig. 3 and Fig. 6 and the `BH_P` time series become visibly less faithful. Since
the conclusion is intensity-dependent, reporting it at one or two intensities is
weaker than reporting the whole trade-off curve (§5.1).

### 4.6 Numerical resolution

`dx = 0.1` in log body mass. The production integral of Eq. (2.4), with its
boundary terms, is only first-order accurate: in this implementation the
discretisation error is about 10% at `dx = 0.1`, falling to 2.5% at `dx = 0.025`
(`tests/test_model.R`). Under `BH_P`, `F_i` is *proportional* to that quantity,
and the error is not the same for every species, so it acts as a species-dependent
bias in the harvest rule. Worth a grid-convergence check on the results, not just
on the rates.

## 5. Three changes that dissolve arbitrariness rather than testing it

### 5.1 Report the trade-off frontier, not a single calibrated point

Instead of picking one `c_P` and comparing, sweep the intensity constant for each
rule over two orders of magnitude and plot total yield against a biodiversity
metric. Each rule becomes a curve. If `BH_P`'s curve dominates the other two over
the whole range, the conclusion is robust to *every* choice of yield-matching
criterion at once, and §3.1 evaporates. If the curves cross, the conclusion holds
only in some yield range, and it matters which.

This is strictly more informative than the paper's design and costs only about
eight runs per rule.

### 5.2 Compare against an unfished control, not against year 0

Removes the confound of §3.3 at negligible cost.

### 5.3 Replicate across assemblages

The paper has one ecosystem for Figs 2–4 and three for Figs 5–6. Between-assemblage
variation in assembly-based models is typically large, and with `n = 4` a
consistent ranking is suggestive, not established. Ten to twenty independent
assemblages would let the frontier be reported with a spread.

## 6. Probably safe

Kernel shape at fixed mean predator:prey mass ratio; the diffusion term; egg mass;
`eps_R`; the exact number of species; `rho_L` and `x_L`; the immigration rate
`I_0` (it only matters where plankton are grazed to near zero). These are worth a
single one-at-a-time check late on, not a campaign.

## 7. The plan

Metrics, computed at year 50 for every run, always relative to the unfished
control rather than to year 0:

* `min_i (B_i/B_i^ctrl)` — the worst-hit species, the paper's implicit measure;
* number of species below 10% and below 1% of control;
* RMS of `log(B_i/B_i^ctrl)` — overall distortion of the assemblage;
* change in Shannon diversity of the biomass vector;
* number of species below the assembly extinction threshold;
* total yield, and cumulative yield over the 50 years.

### Phase 1 — dissolve the calibration question (highest value, lowest cost)

1. **Yield–biodiversity frontier.** For each of the three rules, eight intensity
   multipliers spanning 0.1x to 10x the calibrated constant; 50-year runs; plot
   yield against each metric. Answers §3.1 and §4.5 together.
2. **Unfished control.** Already implemented; fold into every comparison (§3.3).
3. **Alternative calibration criteria.** Recalibrate on year-0 yield and on
   cumulative yield, and report how much `c_P` moves.

*Cost:* about 24 runs per ecosystem at `dt = 0.01`, roughly half an hour each.

### Phase 2 — is the mechanism what is claimed?

4. **Frozen-allocation controls.** `F_i = c_P P_i(0)` held constant, and
   `F_i = c_P/B P_i(0)/B_i(0)` held constant. Separates allocation from feedback
   (§3.5).
5. **`F_i = c_B B_i`.** Tests whether production is doing any work (§3.4).
6. **The `theta` sweep.** `F_i = c B_i^theta` for `theta` in -0.5 to 1.5,
   calibrated to a common yield. Places all rules on one axis.
7. **Predation-release decomposition.** For each rare species, track fishing
   mortality and predation mortality separately over the 50 years, and report how
   much of the biomass change under `BH_P` is attributable to each.

*Cost:* about 20 runs per ecosystem.

### Phase 3 — model assumptions

8. **Reproduction.** `reproduction_level` in {0, 0.1, 0.25, 0.5, 0.75}, refitting
   the steady state each time, then the full three-rule comparison (§3.2). This
   is the most likely single result-changer after §3.1.
8b. **Resource level.** Sweep `I_0/a_0` in {1, 0.1, 0.01, 0} per year and `r_0`
    in {1, 3, 10, 30} g^(1-rho)/yr, reassembling each time; measure the resource
    level across plankton size and the growth elasticity at egg size for each,
    then run the three-rule comparison on the variants where larval competition
    is actually operating (§3.6). Pair with item 8: together these two items
    determine how much recruitment compensation the model has, by either route.

9. **Eq. (A.7).** Repeat the headline comparison under `mu_b_form = "product"`,
   and sweep `mu_b0` and `xi`. Reports how `alpha` in `B ~ P^alpha` moves, since
   the Fig. 4 argument is stated in terms of `alpha` (see `docs/mu_b.md`).
10. **`w_f`.** 100, 200, 400, 800 g; plus a species-relative entry size
    `w_f = w_mat` and `w_f = w_max/4`, which removes the systematic penalty on
    large species (§4.1).
11. **Measurement range for `P` and `B`.** Harvested range (as published) versus
    whole life cycle (§4.2), with `w_f` held at 400 g.
12. **`theta`.** Diagonal:off-diagonal ratios of 1:1, 2.5:1 (published), 5:1;
    plus one heterogeneous matrix (§4.3).
13. **Plankton cap.** 1 g (published) versus 10 g (§4.4).
14. **Grid.** `dx` of 0.1, 0.05, 0.025 on the headline comparison (§4.6).

*Cost:* each variant needs its own assembly, since changing the model changes
which species coexist — roughly 15 minutes of assembly plus half an hour of runs
per variant.

### Phase 4 — how the rule would actually be implemented

15. **Update interval.** `F_i` recomputed every 1, 3, 5 years from the state at
    that moment rather than continuously. Continuous perfect feedback is a strong
    idealisation for a policy claim.
16. **Observation error.** `F_i = c_P * P_i * exp(eta_i(t))` with log-normal
    noise, and with a persistent per-species bias. The paper tests a *fixed*
    per-species intensity factor `z' ~ U(0.5, 1.5)`, which is bias without noise
    and without dynamics; this is the harder test of the robustness-to-uncertainty
    claim.

### Phase 5 — replication

17. Repeat Phases 1 and 2 across 10–20 independently assembled ecosystems and
    report distributions rather than exemplars (§5.3).

### Ordering

Phase 1 first: if the frontier shows `BH_P` dominating everywhere, §3.1
disappears and the rest of the campaign is about *why* rather than *whether*. If
the frontiers cross, everything downstream has to be conditioned on yield level.
Phase 2 next, because it is cheap and tests the stated mechanism directly.
Phase 3 item 8 (recruitment) is the most likely to overturn the headline and
should not be left to the end.

## 8. What each outcome would mean

| Finding | Reading |
|---|---|
| `BH_P` frontier dominates at all yields, in most assemblages | The conclusion is robust and stronger than the paper shows |
| Frontiers cross | The conclusion is a statement about low-intensity fishing; say so |
| Frozen `BH_P` ~ adaptive `BH_P` | The benefit is allocation, not feedback; the policy claim needs rewording |
| `F ∝ B` ~ `F ∝ P` | Production is a proxy; use biomass, which is far easier to estimate |
| Gap closes as `reproduction_level` rises | The result is contingent on the absence of compensatory recruitment |
| Gap closes with species-relative `w_f` | Much of the effect came from a fixed `w_f` penalising large species |
| Most of the rare-species gain is predation release | The mechanism is trophic, not the harvest rule's feedback |
| Gap closes once larval competition is switched on (lower `I_0/a_0` or `r_0`) | The result depends on the recruitment brake being absent, and the paper's own stated coexistence mechanism was not operating at its published parameter values |

## 9. Code hooks needed

Already parameterised: `mu_b_form`, `theta_ii`/`theta_ij` (`LP_FISH`), the
plankton cap (`LP_PLANKTON$w_max`), `dx` (`LP_NUMERICS`), `w_f`
(`lp_set_fishing`), the per-species intensity factor `zf`.

To add:

* new rules in `lp_F()`: `"BHB"` (`c B_i`), `"power"` (`c B_i^theta`), and
  `"frozen"` (a precomputed constant vector of `F_i`);
* a measurement range for `P` and `B` independent of the harvested range — a
  second index alongside `other_params$lp_jf`;
* an update interval and an observation-error model in `lpFMort()` — needs a
  small state object in `other_params`, since mizer's rate functions cannot
  write back to `params`; simplest is to project in annual segments and reset
  the frozen `F_i` between them;
* `reproduction_level(params) <- ...` in `lp_params()`, plus a re-relaxation,
  since imposing `R_max` moves the steady state;
* species-relative `w_f`, which makes the selectivity a full species-by-size
  array rather than one shared mask;
* a mortality-decomposition recorder for item 7, tracking `f_mort` and
  `pred_mort` separately per species over time;
* a resource-level diagnostic — `n_pp / resource_capacity` by size, the growth
  elasticity to fish abundance, and the split of egg-to-recruit survival into
  its growth-mediated and predation-mediated parts (§3.6). The plankton
  parameters themselves are already exposed in `LP_PLANKTON`.
