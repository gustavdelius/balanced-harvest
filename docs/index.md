# Fishing for biodiversity by balanced harvesting

A reimplementation in [mizer](https://sizespectrum.org/mizer/) of

> Law, R. & Plank, M.J. (2023) *Fishing for biodiversity by balanced harvesting.*
> **Fish and Fisheries** 24, 1–17.
> doi:[10.1111/faf.12705](https://doi.org/10.1111/faf.12705)
> (preprint: doi:[10.1101/2021.06.27.450047](https://doi.org/10.1101/2021.06.27.450047))

Source: [github.com/gustavdelius/balanced-harvest](https://github.com/gustavdelius/balanced-harvest)

---

## The question

Can fishing, if scaled sensibly across species, help maintain biodiversity
rather than erode it? Law & Plank compare three ways of setting the fishing
mortality rate of species *i*:

| | rule | |
|---|---|---|
| fixed | $$F_i(t) = F$$ | Eq. (2.7) |
| BH<sub>P</sub> | $$F_i(t) = c_P\,P_i(t)$$ | Eq. (2.8) |
| BH<sub>P/B</sub> | $$F_i(t) = c_{P/B}\,P_i(t)/B_i(t)$$ | Eq. (2.9) |

Every fish enters a single mixed-species fishery at $$w_f = 400\,\mathrm{g}$$,
and everything above that is caught at the same rate, so $$F$$ is
size-independent and the balancing is **across species, not across sizes**.
$$B_i$$ and $$P_i$$ are measured over that same harvested range, and
$$c_{P/B}$$ is simply a fixed exploitation ratio $$E = Y_i/P_i$$ applied to
every species.

The argument is that because biomass and production rate are tightly coupled
across species,

$$
B_i \sim P_i^{\alpha}, \qquad \alpha \approx 1,
$$

setting $$F_i \propto P_i$$ makes fishing mortality roughly proportional to
biomass, so it falls away as a stock falls. Fixed $$F$$ has no such feedback,
and $$P_i/B_i$$ barely changes as a stock thins, so BH<sub>P/B</sub> has almost
none either.

$$\alpha$$ is the regression slope of $$\log B_i$$ on $$\log P_i$$ across the
species of an assemblage, and it recurs throughout this site: it is what makes
BH<sub>P</sub> a feedback rule rather than a fixed allocation, and Figure 4
below is a direct test of it.

## The model, in equations

The quickest way to say what Law & Plank's model is, is to write it down.
Everything here is Appendix A of the paper, restated in mizer's variables: the
per-mass number density $$N_i(w,t)$$ of species *i*, where the paper uses the
log-density $$u_i(x,t) = w\,N_i$$ with $$x = \ln(w/1\,\mathrm{g})$$. The two are
the same PDE divided through by $$w$$; [mapping.md](mapping.md) does that term
by term.

**The dynamics** (Eq. A.1) are McKendrick–von Foerster with a diffusion term:

$$
\frac{\partial N_i}{\partial t}
= -\frac{\partial}{\partial w}\Big[ g_i(w,t)\, N_i \Big]
\;+\; \frac{1}{2}\frac{\partial^2}{\partial w^2}\Big[ D_i(w,t)\, N_i \Big]
\;-\; \mu_i(w,t)\, N_i ,
$$

with eggs entering at $$w_0 = 10^{-3}$$ g at the rate $$R_i(t)$$ given below.
The diffusion term is there because eating a prey item is a discrete jump in
mass: $$D_i$$ is the second moment of the same prey integral whose first moment
gives $$g_i$$.

**Feeding** (Eqs A.2, A.4). Search rate is $$A_i w^{\alpha_q}$$, and the feeding
kernel is a **box** spanning four decades of predator–prey mass ratio,
$$e^{\beta \pm 3\sigma} = 10$$ to $$10^5$$, with height $$1/(6\sigma)$$:

$$
E_i(w,t) = \frac{A_i\, w^{\alpha_q}}{6\sigma}
  \int_{w/10^{5}}^{\,w/10}
  \Big[\, N_R(w_p,t) + \sum_j \theta_{ij} N_j(w_p,t) \Big]\, w_p\,
  \mathrm{d}w_p .
$$

There is **no satiation and no metabolic cost**, so the whole assimilated
intake is split between growth and reproduction by the allocation
$$\epsilon_i(w) = 1 - \psi_i(w)$$ of Eq. (A.8):

$$
g_i(w,t) = \epsilon_i(w)\, K\, E_i(w,t),
\qquad
\psi_i(w) = \frac{1}{1 + (w/w_{\mathrm{mat},i})^{-\rho_m}}
            \left(\frac{w}{w_{\max,i}}\right)^{\rho_\infty} ,
$$

with $$w_{\mathrm{mat},i} = w_{\max,i}/10$$. Eggs are produced without any
stock-recruitment relationship (Eq. A.9):

$$
R_i(t) = \frac{\epsilon_R}{w_0}\int \psi_i(w)\, K\, E_i(w,t)\, N_i(w,t)\,
         \mathrm{d}w .
$$

**Mortality** is predation, a larval term and a background term,
$$\mu_i = \mu^{\mathrm{pred}}_i + \mu^{L}_i + \mu^{b}_i$$. Predation
(Eq. A.3) is the mirror of the encounter integral,

$$
\mu^{\mathrm{pred}}_i(w_p,t) = \sum_j \theta_{ji}\,
  \frac{A_j}{6\sigma} \int_{10\,w_p}^{10^{5} w_p}
  w^{\alpha_q}\, N_j(w,t)\, \mathrm{d}w ,
$$

the larval term (Eq. A.6) is a reverse sigmoid that has died away by
$$w_L = 0.1$$ g,

$$
\mu^{L}_i(w) = \frac{\mu_{\mathrm{egg},i}}{1 + (w/w_L)^{\rho_L}} ,
$$

and the background term (Eq. A.7) is **food-dependent** — metabolic need
relative to how fast food is arriving, so it is a starvation feedback rather
than a constant:

$$
\mu^{b}_i(w,t) = \mu_b^{(0)} \left(\frac{w}{w_0}\right)^{-\xi}
  \frac{\tilde g_i(w_0,t)}{\tilde g_i(w,t)},
  \qquad \tilde g_i = K E_i / w .
$$

That last ratio is printed the other way up in the paper, which changes the
model qualitatively; [mu_b.md](mu_b.md) sets out why this implementation
inverts it.

**The plankton** (Eqs A.11-A.14) is logistic with immigration, and has no size
structure of its own — cells neither grow nor are born into a smaller class:

$$
\frac{\partial N_R}{\partial t}
= r(w)\, N_R \left(1 - \frac{N_R}{a(w)}\right)
\;-\; \mu^{\mathrm{pred}}_R(w,t)\, N_R + I(w),
$$

with $$r(w) = r_0 w^{-\rho}$$, and carrying capacity $$a$$ and immigration
$$I$$ both power laws in $$w$$ anchored at 1 mg with $$I_0 = a_0$$ — so
$$I/a = 1$$ yr⁻¹ at *every* plankton size, which turns out to matter a great
deal — see [Is the result robust?](#is-the-result-robust) below.

**What is harvested** (Eqs 2.3, 2.4, 2.6). Everything above the entry mass
$$w_f$$ is fished at the same rate, and biomass and production are measured
over that same range:

$$
B_i(t) = \int_{w_f}^{w_{\max,i}} w\, N_i(w,t)\, \mathrm{d}w,
\qquad
Y_i(t) = F_i(t)\, B_i(t),
$$

$$
P_i(t) = \underbrace{\int_{w_f}^{w_{\max,i}} g_i(w,t)\, N_i(w,t)\,
\mathrm{d}w}_{\text{growth inside the range}}
\;+\; \underbrace{\big[w\, g_i N_i\big]_{w_f}}_{\text{flux in}}
\;-\; \underbrace{\big[w\, g_i N_i\big]_{w_{\max,i}}}_{\text{flux out}} .
$$

$$P_i$$ is **somatic production** — the rate at which the harvested part of the
stock puts on flesh, plus the biomass that grows across $$w_f$$ into the range.
It is not yield, and not surplus production in the stock-assessment sense. The
flux out of the top vanishes, because all food goes to reproduction at
$$w_{\max}$$.

**The parameters** are the paper's Tables 1 and 2, reproduced verbatim in
[`R/lp_constants.R`](https://github.com/gustavdelius/balanced-harvest/blob/main/R/lp_constants.R):

| | | |
|---|---|---|
| $$A_i$$ | 37.5 m² yr⁻¹ g$$^{-\alpha_q}$$ | search-rate coefficient |
| $$\alpha_q$$ | 0.85 | search-rate exponent (mizer's `q`) |
| $$\beta,\ \sigma$$ | 6.908, 1.535 | centre and width of the box kernel, in log PPMR |
| $$K$$ | 0.1 | conversion efficiency (mizer's `alpha`) |
| $$\mu_b^{(0)},\ \xi$$ | 0.1 yr⁻¹, 0.15 | background mortality at egg size, and its size scaling |
| $$w_L,\ \rho_L$$ | 0.1 g, 5 | where larval mortality dies away, and how sharply |
| $$\rho_m,\ \rho_\infty$$ | 15, 0.2 | sharpness of maturation, approach to $$w_{\max}$$ |
| $$\epsilon_R$$ | 0.2 | reproductive efficiency |
| $$\theta$$ | 0.5 / 0.2 / 1 | cannibalism / between species / on plankton (Appendix C) |
| $$r_0,\ \rho$$ | 10 yr⁻¹ at 1 g, 0.15 | plankton rate of increase and its size scaling |
| $$a_0 = I_0$$ | 2000 m⁻² at 1 mg | plankton carrying capacity and immigration |
| $$\lambda$$ | 2 | slope of the isolated plankton spectrum |
| $$w_f$$ | 400 g | entry mass of the fishery |
| $$\Delta x,\ \Delta t$$ | 0.1, 0.002 yr | log-mass grid step, time step |

Species differ in exactly two things: maximum body mass $$w_{\max,i}$$ and
larval death rate $$\mu_{\mathrm{egg},i}$$. Everything else is shared across
species. That single fact is what [life-history.md](life-history.md) is about.

## Is this really the paper's model?

The paper's model is **not** a standard mizer model. It has no satiation and no
metabolic cost, no imposed stock-recruitment relationship, food-dependent
intrinsic mortality, logistic plankton with immigration instead of a
semichemostat resource, and a box feeding kernel spanning four decades of prey
size. All of that is built from Appendices A–C; the term-by-term correspondence
is in [mapping.md](mapping.md).

Rather than assert the mapping, [`tests/test_model.R`](https://github.com/gustavdelius/balanced-harvest/blob/main/tests/test_model.R)
checks mizer's rates against hand-coded integrations of the paper's own
equations on the model's own grid. All sixteen checks pass — the encounter rate
against Eq. (A.2) and predation mortality against Eq. (A.3) agree to about
1e-13, and the production integral of Eq. (2.4), including its boundary terms,
converges at the expected first order.

Three independent quantities then came out close to the paper's without being
fitted to them:

| | this reimplementation | Law & Plank |
|---|---|---|
| $$\alpha$$ in $$B \sim P^{\alpha}$$ (four ecosystems) | **0.977, 0.998, 0.993, 0.984** | 1.004 |
| $$c_P$$ calibrated against $$F = 0.1$$ | **0.906** m² g⁻¹ | 1 |
| $$c_{P/B}$$ calibrated against $$F = 0.1$$ | **0.259** | 0.25 |

$$\alpha$$ matters most: it is the exponent the whole Figure 4 argument rests
on ("species were near a line of slope $$1+\alpha = 2.004$$"), it is set by how much
adult mortality the model has, and it emerges here from four independently
assembled ecosystems.

## The ecosystems

There is no species pool taken from nature here. Following Appendix B, each
ecosystem is **assembled inside the model**: introduce one fish species at a
time, with maximum body mass drawn log-uniform over 100 g – 40 kg and larval
death rate $$\mu_\mathrm{egg}$$ drawn uniform over 28–32 yr⁻¹, start it at a
low density on a power law, let the system relax for 50 years, cull anything
below 2×10⁻⁶ g m⁻², and repeat until 15 species coexist. An assemblage is
accepted only if it meets Appendix B's four criteria: the species coexist, the
state is close to equilibrium, $$w_{\max}$$ spans a wide range inside the drawn
limits, and biomass spans roughly four orders of magnitude.
[`tests/test_ecosystem.R`](https://github.com/gustavdelius/balanced-harvest/blob/main/tests/test_ecosystem.R)
checks all four, plus the numbers the paper quotes in Sections 2.4 and 3.

Four ecosystems are used throughout this site, and their names recur on every
page:

- **`eco1`** — seed 101, every species given the same search-rate coefficient
  $$A_i$$. This is the one Figures 2, 3 and 4 are drawn from, and the one most
  of the experiments run on.
- **`eco_r1`, `eco_r2`, `eco_r3`** — seeds 201–203, with $$A_i$$ randomised per
  species ($$A_i = z_i A$$, $$z_i \sim N(1, 0.1)$$), which is what the paper does
  for its Figures 5 and 6. The extra life-history variation is the only
  difference.

Later pages add more, all named after the seed they were assembled from:
`eco301`–`eco312` are twelve further ecosystems assembled under exactly the
`eco1` protocol, used for the replication in
[robustness-results.md](robustness-results.md); and the `varied_*`,
`fixed_*`, `neutral_*` and `gradient_*` runs of
[life-history.md](life-history.md) are named *experimental arm* underscore
*seed*.

| | eco1 | eco_r1 | eco_r2 | eco_r3 | paper |
|---|---|---|---|---|---|
| species | 15 | 15 | 15 | 15 | 15 |
| invasion attempts | 23 | 25 | 30 | 33 | 25, 24, 19, 45 |
| $$w_{\max}$$ range | 327 g – 13.6 kg | 860 g – 9.3 kg | 459 g – 14.5 kg | 566 g – 22.3 kg | 100 g – 40 kg |
| biomass span | 5.8 decades | 3.0 | 3.8 | 7.8 | ~4 |
| cor(log B, log P) | 0.999 | 0.998 | 0.997 | 1.000 | — |
| primary production (g m⁻² yr⁻¹) | 4310 | 4332 | 4324 | 4312 | ~4000 |
| total fish biomass (g m⁻²) | 5.95 | 6.24 | 5.82 | 5.51 | ~4–5 |
| yield at `F = 0.1` (g m⁻² yr⁻¹) | 0.287 | 0.503 | 0.416 | 0.289 | ~0.25 |
| Fogarty ratio (‰) | 0.067 | 0.116 | 0.096 | 0.067 | 0.06 |

The paper does not report its random seed or the realised life histories of its
15 species, so this reproduces the assembly *procedure*, not the particular
assemblage. The figures below are therefore comparable in kind, not number for
number.

---

## Figure 2 — the biomass–production relationship

![Figure 2](figures/fig2.png)

Panel (b) is the relationship the whole argument depends on: biomass and
production rate are tightly correlated across species, over six decades. Panel
(c) shows the pattern the paper reports in its own Fig. 2c — the small-$$w_{\max}$$
species are the common ones, and the rare species sit at larger maximum body
mass.

Why panel (b) comes out this way, and what it takes to break it, is worked out
in [life-history.md](life-history.md): $$P/B$$ is the biomass-weighted mean
mortality rate, and mortality here is a function of body size that no species
can escape, so abundance cancels exactly. Two attempts to spread $$P/B$$ across
species are reported there. Giving every species its own egg mass fails —
assembly drives all survivors to the largest egg on offer. Giving them their
own activity succeeds, but only once faster species are made correspondingly
more vulnerable to predation, which turns activity into a neutral axis that
assembly cannot climb; a genuine fast-slow continuum then assembles, with a
ten-fold range in age at maturity. $$B \sim P$$ survives even that.

Running the three rules on those ecosystems is also the sharpest test of the
paper's conclusion so far. BH<sub>P/B</sub> is no longer a near-duplicate of a
fixed $$F$$ there — it spans an order of magnitude and tracks turnover at
$$r = 0.94$$, which is orthodox single-stock advice — and it *still* leaves the
rare species exposed, because turnover is orthogonal to rarity. Only tracking
abundance protects biodiversity.

(The paper's Fig. 2a is a scatter from an Ecopath model of the West Scotland
shelf. That is empirical data, not model output, and is not reproduced here.)

## Figure 3 — three ways to harvest the same assemblage

![Figure 3](figures/fig3.png)

All three regimes are calibrated to the same total yield after 50 years
(0.197–0.198 g m⁻² yr⁻¹).

**Top row**, total biomass of each species over the 50 years of fishing. Under
fixed `F` and BH<sub>P/B</sub> several species fall along straight lines in log
biomass — exponential decline. Under BH<sub>P</sub> the trajectories are
comparatively flat.

**Middle row** is the mechanism, laid bare. Fixed `F` is flat across species by
construction. BH<sub>P/B</sub> is *also nearly flat* — which is the paper's
point: `P_i/B_i` is much the same for every species, so a constant exploitation
ratio is barely distinguishable from a constant `F`. BH<sub>P</sub> spans seven
orders of magnitude, dropping to `F ≈ 10⁻⁷ yr⁻¹` for the rarest species.

**Bottom row**, the ratio of year-50 to year-0 biomass. Under fixed `F` and
BH<sub>P/B</sub> the larger-numbered (rarer) species fall well below the dotted
line; under BH<sub>P</sub> they stay close to it.

## Figure 4 — yield against production rate

![Figure 4](figures/fig4.png)

Dotted lines are isoclines of constant exploitation ratio `E`; the dashed line
in the middle panel has slope $$1+\alpha$$.

Under fixed `F` and BH<sub>P/B</sub> the species lie along a line of slope 1 —
every species is exploited at the same ratio, regardless of how rare it is. That
is what puts the rare ones on a downward path. Under BH<sub>P</sub> the species
lie along the steeper slope-$$(1+\alpha)$$ line: yield rises faster than linearly with
production, so rare species are exploited at a progressively *smaller* fraction
of their production. The open circle and arrow trace the rarest fished species
from year 0 to year 50.

This is the paper's central quantitative claim, and it reproduces.

## Figure 5 — emergent growth

![Figure 5](figures/fig5.png)

Growth trajectories obtained by integrating Eq. (C.1). Nothing here is imposed:
body growth emerges from what the fish actually eat. By age 1 the species span a
3.3- to 3.7-fold range in body mass, against the 4- to 7-fold range the paper
reports for its own assemblages.

## Figure 6 — three more ecosystems, harder fishing

![Figure 6](figures/fig6.png)

Three independently assembled ecosystems (rows) with extra life-history
variation, baseline fishing doubled to 0.2 yr⁻¹, and a random per-species
intensity factor $$z'_i \sim U(0.5, 1.5)$$ applied identically to all three rules. The
BH<sub>P</sub> column is visibly flatter than either neighbour in all three
assemblages, and the BH<sub>P/B</sub> collapses in assemblage 3 are dramatic.

---

## The headline comparison

On `eco1`, with all three regimes calibrated to the same total yield:

| regime | yield (g m⁻² yr⁻¹) | worst species, vs year 0 | worst species, vs unfished control | species below 10% of control |
|---|---|---|---|---|
| unfished control | 0 | 0.174 | 1.00 | 0 |
| fixed *F* | 0.198 | 0.046 | 0.166 | 0 |
| **BH<sub>P</sub>** | **0.198** | **0.345** | **0.573** | **0** |
| BH<sub>P/B</sub> | 0.197 | 0.025 | 0.070 | 3 |

A caveat on the "worst species" columns, found while running these rules on the
fast-slow ecosystems of [life-history.md](life-history.md): under
BH<sub>P</sub> the worst-affected species is usually the *most abundant* one,
because that is where the rule deliberately concentrates effort. The columns
therefore measure harm to a common species under one rule and harm to a rare
one under another. Conditioning on rarity instead changes BH<sub>P</sub>'s
measured advantage there from 0.95-1.6x to 5.4-16.9x.

Across all four ecosystems, counting species that end below 10% and below 1% of
what they would have been with no fishing at all:

| | eco1 | eco_r1 | eco_r2 | eco_r3 |
|---|---|---|---|---|
| fixed *F* | 0 / 0 | 0 / 0 | 3 / 0 | 5 / 2 |
| **BH<sub>P</sub>** | **0 / 0** | **0 / 0** | **0 / 0** | **0 / 0** |
| BH<sub>P/B</sub> | 3 / 0 | 3 / 1 | 6 / 3 | 6 / 6 |

BH<sub>P</sub> does not push a single species below 10% of its unfished
trajectory in any of the four ecosystems, and its worst-affected species retains
25–57% of unfished biomass. In the three Figure 6 ecosystems it does this while
returning *more* yield than fixed `F` (0.797 vs 0.631, 0.539 vs 0.422, 0.369 vs
0.344), so it is not buying protection by fishing less.

**The paper's conclusion reproduces.**

## One correction to how the comparison is framed

Law & Plank compare year-50 fished biomass against year-0 unfished biomass. But
the assembled ecosystem is explicitly only a quasi-equilibrium — the paper says
so — and it drifts. In the unfished control run here, the worst-affected species
still falls to **17% of its year-0 biomass over 50 years with no fishing at
all**, and individual species move by 39–86% over a further 20 unfished years.

Comparing against year 0 therefore charges fishing for drift it did not cause.
The effect is not small: fixed `F` on `eco1` looks like it takes its worst
species to 4.6% of where it started, but only to 17% of where it would have been
anyway, and the count of species below 10% falls from 2 to 0 once the control is
used. The ranking of the three rules is unchanged, but the magnitudes are not.
Both columns are reported above.

## Where this departs from the paper

Two places in the paper needed a judgement call, and one of them matters a great
deal. Full detail in [mapping.md](mapping.md) and [mu_b.md](mu_b.md).

**Eq. (A.7) contradicts itself.** The prose says background mortality is
proportional to metabolic need *divided by* the rate food becomes available; the
printed formula *multiplies* by it. The two differ by two orders of magnitude in
adult mortality. The printed version gives near-zero adult mortality, biomass
piling up at $$w_{\max}$$, roughly ten times too much fish biomass, and large species
so dominant that small ones cannot invade at all. The prose version reproduces
the paper's own reported biomass range, its 0.25 g m⁻² yr⁻¹ yield, its 0.06 ‰
Fogarty ratio, and its ordering in which small species are common and large ones
rare. This implementation follows the prose; `mu_b_form = "product"` gives the
other.

Three smaller ones: Eq. (A.3) as printed evaluates the predator's search area at
the *prey's* body mass (Eq. A.2 uses the consumer's own size, so this is a typo);
Table 1 contradicts the text on the plankton carrying-capacity offset (the text
is right, and gives the stated 4000 g m⁻² yr⁻¹ of primary production); and
Appendix B does not say what power law an invading spectrum starts on.

## Is the result robust?

Reproducing a result is not the same as testing it.
[robustness.md](robustness.md) works through which of the paper's arbitrary
choices could change the conclusion — the year-50 yield-matching criterion, the
absence of compensatory recruitment, whether `F ∝ B` would do just as well as
`F ∝ P`, and whether the *adaptive feedback* or merely the *initial allocation*
is doing the work — and sets out an investigation.

**[robustness-results.md](robustness-results.md) reports what came out of it** —
phases 1, 2, 4 and 5 and Phase 3 items 8, 8b, 10 and 11, across up to thirteen
independently assembled ecosystems.
In short: the calibration criterion turns out not to matter (the frontiers never
cross, and I was wrong to call it the most load-bearing choice); $$F$$
proportional to *biomass* does everything $$F$$ proportional to production does;
the adaptive feedback turns out to be essential, though it took thirteen
ecosystems to establish that — on one it looked as though the initial allocation
did most of the work; the plankton cannot be made contestable
without destroying the large-bodied community the paper harvests; and, most
consequentially, imposing a stock-recruitment relationship on the published
ecosystem while holding the community exactly fixed shrinks
BH<sub>P</sub>'s advantage over fixed $$F$$ from 3.4× to 1.0×. Replicated over
thirteen independently assembled ecosystems, the ranking
BH<sub>P</sub> > fixed > BH<sub>P/B</sub> holds 13 times out of 13, with the
advantage varying fourfold (median 2.4×). The feedback tolerates being updated
only once a decade, but not observation error much beyond a factor of two.

One finding from that work is already in hand. Appendix B nominates competition
among larvae for plankton as the mechanism that replaces an imposed
stock-recruitment relationship. Measured on the assembled 15-species ecosystem,
the plankton over the entire larval prey range (10⁻⁸ to 10⁻⁴ g) sits **at its
carrying capacity**, and doubling the whole fish spectrum changes growth at egg
size by only −1.3% (against −91% at 100 g). The reason is structural: the paper
sets immigration `I₀` equal to carrying capacity `a₀` with the same size
scaling, so `I/a = 1 per year at every plankton size`, while plankton at larval
prey sizes regenerate at 40–220 yr⁻¹. The model therefore has neither an imposed
recruitment brake nor, at these parameter values, the emergent one offered in
its place.

## Reproducing this

```bash
Rscript tests/test_model.R        # verify the mapping; needs no inputs
Rscript run_assembly.R            # ~1 h  -> data/ecosystems.rds
Rscript tests/test_ecosystem.R    # check the ecosystems against Appendix B
Rscript run_figures.R             # ~40 m -> data/results.rds, figures/
Rscript tests/test_convergence.R  # the numerical liberties taken
```

The life-history experiments of [life-history.md](life-history.md) are separate
and slower, one assembly per invocation:

```bash
Rscript run_egg_assembly.R 301 varied        # -> data/egg/
Rscript run_activity_assembly.R 301 neutral  # -> data/activity/
Rscript run_activity_harvest.R 301           # -> data/activity_harvest/
```

Requires [mizer](https://sizespectrum.org/mizer/) (developed against 3.3.0.9000),
ggplot2 and patchwork.

## Where does B ~ P come from?

[life-history.md](life-history.md) works out why production is proportional to
biomass in these ecosystems, and tries twice to break the relationship by
giving species life histories that vary independently of maximum body mass.

## A second test, on a real ecosystem

The ecosystems above are assembled inside the model, which is the right way to
test the paper but says nothing about whether the result survives contact with a
real one. [north-sea.md](north-sea.md) runs the same three rules on mizer's
North Sea model with real depleted species — thornback ray and Atlantic halibut
— added at abundances estimated from ICES IBTS survey data.

The conclusion holds there too, but conditionally: `NS_params` as it ships has
nothing capable of being lost, and BH<sub>P</sub> looks *worst* of the three
until something at risk is put in. Halibut then goes from 2% of its unfished
trajectory under a constant `F` to 540% under BH<sub>P</sub>, at equal yield,
while BH<sub>P/B</sub> destroys the thornback ray.

## Caveats

The ecosystems are not the paper's ecosystems — the assembly procedure is
reproduced, not the particular assemblage, because the seed and realised life
histories are not reported. Four ecosystems is a small sample for a claim about
assemblages in general. And the comparison rests on a single calibration
criterion; the trade-off frontier proposed in [robustness.md](robustness.md)
would be a stronger design than the one used here and in the paper.

