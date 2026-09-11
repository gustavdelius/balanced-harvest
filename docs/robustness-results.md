# Robustness investigation: results

What came out of the plan in [robustness.md](robustness.md). Phases 1 and 2 and
Phase 3 items 8 and 8b are complete. Everything below is measured, not argued.

Three headline results:

1. The calibration criterion the plan called "the single most load-bearing
   arbitrary choice" does not matter at all. The rules' yield-biodiversity
   frontiers never cross.
2. Production is not the active ingredient. $$F$$ proportional to *biomass* does everything $$F$$ proportional to
   production does, and about three-quarters of
   the benefit comes from the initial allocation across species rather than
   from the adaptive feedback the paper emphasises.
3. The plankton cannot be made contestable without destroying the large-bodied
   community the paper harvests - and in the communities that survive doing so,
   BH<sub>P</sub>'s advantage inverts.
4. Imposing a stock-recruitment relationship on the published ecosystem, holding
   the community exactly fixed, shrinks BH<sub>P</sub>'s advantage from 3.4x to
   1.0x. The paper's result is contingent on its deliberate choice to have no
   compensatory recruitment.

Reproduce with `Rscript run_robustness.R` and `Rscript run_resource_sweep.R`.

---

## The experiments, stated precisely

Notation follows the paper: species $$i = 1,\dots,n$$, body mass $$w$$, number
density $$N_i(w,t)$$ per unit mass per m². The model grid is
$$w_1 < \dots < w_J$$ with $$w_1 = 1$$ mg and bin widths $$\Delta w_j$$, and
$$g_{ij}(t)$$ is mizer's `e_growth`, the somatic growth rate $$dw/dt$$ in
g yr⁻¹ — equal to $$\epsilon_i(w)\,\tilde g_i(w,t)$$ in the paper's Appendix A
notation. Every experiment uses $$T = 50$$ years of fishing.

### The harvest rules

The fishery is knife-edged at a shared entry mass $$w_f$$, taken as the first
grid point at or above it,

$$
j_f = \min\{\, j : w_j \ge w_f \,\}.
$$

All rules have the same form,

$$
F_i(w_j, t) = F_i(t)\,\mathbf{1}[\,j \ge j_f\,],
\qquad
F_i(t) = c\, z'_i\, g_i(t),
$$

differing only in the per-species allocation $$g_i(t)$$:

$$
g_i =
\begin{cases}
1 & \text{fixed, Eq. (2.7)}\\[2pt]
P_i(t) & \mathrm{BH}_P,\ \text{Eq. (2.8)}\\[2pt]
P_i(t)/B_i(t) & \mathrm{BH}_{P/B},\ \text{Eq. (2.9)}\\[2pt]
B_i(t) & \mathrm{BH}_B\\[2pt]
B_i(t)^{\theta} & \mathrm{power}(\theta)\\[2pt]
h_i & \mathrm{frozen}(h),\ \text{constant in } t
\end{cases}
$$

Negative or non-finite $$F_i$$ are set to zero. $$z'_i$$ is the per-species
fishing intensity factor: $$z'_i = 1$$ for eco1, and
$$z'_i \sim U(0.5, 1.5)$$ (seed 99) for eco_r1–3, held fixed across rules within
an ecosystem, as in the paper's Fig. 6.

The frozen controls take $$h_i = g_i^{X}(0)$$, the year-0 allocation of rule
$$X$$, and hold it constant for the whole 50 years. Comparing
$$\mathrm{frozen}(g^{\mathrm{BH}_P}(0))$$ with $$\mathrm{BH}_P$$ separates
*allocating* effort across species from *updating* that allocation as the
ecosystem changes.

### $$B_i$$ and $$P_i$$, as discretised

Over the harvested range $$[w_f, w_{\max,i}]$$, with $$J$$ the top grid index:

$$
B_i(t) = \sum_{j \ge j_f} w_j\, N_{ij}(t)\, \Delta w_j
\tag{2.3}
$$

$$
P_i(t) = \underbrace{\sum_{j \ge j_f} g_{ij}(t)\, N_{ij}(t)\, \Delta w_j}_{\text{production inside the range}}
\;+\; \underbrace{w_{j_f}\, g_{i j_f}(t)\, N_{i j_f}(t)}_{\text{flux in at } w_f}
\;-\; \underbrace{w_{J}\, g_{iJ}(t)\, N_{iJ}(t)}_{\text{flux out at the top}}
\tag{2.4}
$$

The boundary terms are the ones integration by parts demands; the flux out
vanishes because $$\epsilon_i \to 0$$ at $$w_{\max,i}$$. Yield is
$$Y_i(t) = F_i(t)\,B_i(t)$$, Eq. (2.6), which is exact because $$B_i$$ uses the
same $$j_f$$ as the selectivity.

$$B_i^{\mathrm{tot}}(t)$$ denotes the same sum over the *whole* grid,
$$j \ge 1$$. The biodiversity measures below use $$B^{\mathrm{tot}}$$; the
yields use $$B_i$$.

### Anchoring the intensity constant

Each rule's constant is set so that the *initial* total yield equals that of
fixed fishing at a reference rate $$F_{\mathrm{ref}}$$. Since
$$Y(0) = \sum_i F_i(0) B_i(0) = c \sum_i z'_i g_i(0) B_i(0)$$,

$$
c^0_X \;=\; F_{\mathrm{ref}}\,
  \frac{\sum_i B_i(0)}{\sum_i z'_i\, g^X_i(0)\, B_i(0)},
$$

which is closed-form and needs no simulation.
$$F_{\mathrm{ref}} = 0.1\,\mathrm{yr}^{-1}$$ for eco1 and the plankton variants,
and $$0.2\,\mathrm{yr}^{-1}$$ for eco_r1–3, matching the paper's doubled
baseline for its Fig. 6. The sweep then runs $$c = m\,c^0_X$$ for

$$
m \in \{0.1,\, 0.25,\, 0.5,\, 1,\, 2,\, 4,\, 8\}
\quad\text{(Phases 1, 2 and the recruitment sweep)},
$$

$$
m \in \{0.25,\, 0.5,\, 1,\, 2,\, 4\}
\quad\text{(plankton variants)}.
$$

### Outcome measures

The unfished control is the same params object with $$F_i \equiv 0$$, started
from the same state and run with the same $$\Delta t$$ and $$T$$, so that
numerical and drift effects cancel. Writing

$$
R_i = \frac{B_i^{\mathrm{tot}}(T)}{B_i^{\mathrm{tot,\,ctrl}}(T)},
$$

the measures are

$$
\begin{aligned}
\text{worst species} &= \min_i R_i, \\
\text{below 10\%} &= \#\{\, i : R_i < 0.1 \,\}, \\
\text{RMS log} &= \sqrt{\tfrac1n \textstyle\sum_i \big(\log \max(R_i, 10^{-12})\big)^2}, \\
\text{terminal yield} &= \textstyle\sum_i F_i(T)\, B_i(T), \\
\text{cumulative yield} &= \textstyle\int_0^T \sum_i Y_i(t)\, dt
   \quad\text{(trapezoid on the saved times)}.
\end{aligned}
$$

The $$10^{-12}$$ floor only binds for species already extinct by any standard.

### Comparing rules at a common yield

Each rule gives a set of (terminal yield, metric) points along its frontier. A
metric is read at a yield none of the runs hit exactly by interpolating linearly
in $$\log(\text{yield})$$; log–log if the metric is strictly positive along the
frontier, otherwise linear in the metric (which is what the integer "below 10%"
count uses). The reference yield is always the fixed-$$F$$ run at $$m = 1$$ in
that ecosystem.

### Plankton variants

In mizer's per-mass densities the paper's Eqs (A.12)–(A.14) read

$$
r(w) = r_0\, w^{-\rho},
\qquad
a(w) = \frac{a_0}{w}\left(\frac{w}{w_a}\right)^{1-\lambda},
\qquad
I(w) = \frac{I_0}{w}\left(\frac{w}{w_I}\right)^{1-\lambda},
$$

with $$\rho = 0.15$$, $$\lambda = 2$$, $$w_a = w_I = 1$$ mg. Total primary
production is $$\int r(w)\, w\, N_{pp}(w)\, dw$$, which for an undepleted
spectrum is proportional to $$r_0 a_0$$. The variants therefore hold

$$
r_0\, a_0 = 2\times 10^{4}
\qquad \text{(the published pair is } 10,\ 2000\text{)},
$$

so that they differ in how *contested* the plankton is, not in how much of it
there is. Variants used, as $$(r_0, a_0, I_0)$$: $$(10,\, 2\!\times\!10^3,\,
2\!\times\!10^3)$$ published, $$(3,\, 6.67\!\times\!10^3,\, 2\!\times\!10^3)$$,
$$(1,\, 2\!\times\!10^4,\, 2\!\times\!10^3)$$, and in the single-species screen
also $$(0.5,\, 4\!\times\!10^4,\, 2\!\times\!10^3)$$ and
$$(0.1,\, 2\!\times\!10^5,\, 0)$$.

### Resource level

Writing $$L = N_{pp}/a$$ and $$\iota = I/a$$, the steady state of Eq. (A.11) at
grazing mortality $$d$$ is

$$
\iota + r\,L\,(1-L) - d\,L = 0
\qquad\Longrightarrow\qquad
r\,L^2 - (r-d)\,L - \iota = 0,
$$

$$
L = \frac{(r-d) + \sqrt{(r-d)^2 + 4r\iota}}{2r}.
$$

This is the same root that `lpResource()` steps towards, so a model at
equilibrium reproduces it exactly.

### Strength of larval competition

The elasticity of somatic growth to fish abundance, at fixed spectrum shape,
letting the plankton re-equilibrate:

$$
\varepsilon(w) = \frac{
  \ln g\big(w;\, 2N,\, N_{pp}^{*}(2N)\big) -
  \ln g\big(w;\, N,\, N_{pp}^{*}(N)\big)}{\ln 2},
$$

where $$N_{pp}^{*}(N)$$ is the root above evaluated with the grazing mortality
$$d$$ that the fish abundance $$N$$ generates. Reported as the mean over species
at $$w = w_1$$ (egg size) and $$w = 1$$ g. It is zero if growth is independent
of how many fish are competing, which is what "larval competition is inert"
means.

### Matching the fishery across variants

Contesting the plankton changes the community's size structure so much that a
common $$w_f$$ in grams is not a common experiment. Instead $$w_f$$ is chosen
per variant so that the same *fraction of community biomass* is exposed:

$$
\varphi(w) = \frac{\sum_{j:\, w_j \ge w} \sum_i w_j N_{ij} \Delta w_j}
                  {\sum_{j} \sum_i w_j N_{ij} \Delta w_j},
\qquad
w_f^{X} = \arg\min_{w} \big|\varphi_X(w) - \varphi_{\text{pub}}(400\,\mathrm{g})\big|,
$$

with $$\varphi_{\text{pub}}(400\,\mathrm{g}) = 0.483$$, giving $$w_f$$ of 402 g,
10 g and 1 g for $$r_0 = 10,\ 3,\ 1$$.

### Compensatory recruitment

`setBevertonHolt()` imposes

$$
R_{dd} = \frac{R_{di}\, R_{\max}}{R_{di} + R_{\max}}
$$

while rescaling `erepro` so that the initial state is preserved exactly. With
the reproduction level $$L = R_{dd}/R_{\max}$$ evaluated at the initial state,

$$
R_{\max} = R_{di}\,\frac{1-L}{L},
\qquad
\frac{d \ln R_{dd}}{d \ln R_{di}} = \frac{R_{\max}}{R_{di} + R_{\max}} = 1 - L,
$$

so $$L = 0$$ is the paper — recruitment tracks egg production one for one — and
$$L \to 1$$ is recruitment pinned at $$R_{\max}$$ regardless of spawning stock.
Preserving the steady state requires
$$\mathrm{erepro}(L) = \mathrm{erepro}(0)/(1-L)$$.

### Numerics

Sweeps and assembly use $$\Delta t = 0.01$$ yr; the final relaxation of each
assembled ecosystem uses the paper's $$\Delta t = 0.002$$ yr.
`tests/test_convergence.R` checks that the two agree to 0.3% and that the
closed-form plankton step agrees with the paper's explicit Euler to 0.02%.
Saved time steps are every 5 years. Stability is judged from an unfished
200-year projection as the coefficient of variation of total biomass over years
100–200; limit-cycle amplitude and period come from a 600-year projection, using
years 300–600 and the mean spacing of successive maxima.

---

## Phase 1 — the calibration criterion does not matter

I had called the paper's choice to match total yield *after 50 years* "the
single most load-bearing arbitrary choice". **That was wrong.**

Recalibrating against cumulative yield over the 50 years instead of terminal
yield moves the constants by about 2%:

| | terminal yield | cumulative yield |
|---|---|---|
| $$c_P$$ | 0.918 | 0.901 |
| $$c_{P/B}$$ | 0.263 | 0.258 |

And the stronger version of the test — sweeping each rule's intensity over an
80-fold range and comparing whole frontiers rather than single points — shows
the frontiers **never cross**. At matched terminal yield, the worst-affected
species retains (as a fraction of its unfished trajectory):

| ecosystem | fixed *F* | BH<sub>P</sub> | BH<sub>P/B</sub> |
|---|---|---|---|
| eco1 | 0.166 | **0.566** | 0.061 |
| eco_r1 | 0.106 | **0.490** | 0.012 |
| eco_r2 | 0.0031 | **0.345** | 7e-6 |
| eco_r3 | 0.0104 | **0.297** | 5e-7 |

The shape matters as much as the ranking. Over an 80-fold increase in fishing
intensity, BH<sub>P</sub>'s worst species declines gently, from 0.85 to 0.09,
while fixed *F* and BH<sub>P/B</sub> fall super-exponentially — in eco_r2, from
0.58 to 1e-18 and from 0.51 to below 1e-100. That is precisely the density
dependence the paper describes: as a stock thins its own fishing mortality
falls, so it cannot be driven to zero, whereas under a constant *F* or a
constant exploitation ratio nothing arrests the decline.

**So the conclusion is robust to the calibration choice, and more robustly true
than the paper's single-point comparison shows.** Section 3.1 of the plan is
answered and closed.

One qualification. On the whole-assemblage distortion measure (RMS of
`log(B/B_control)`), BH<sub>P</sub> is better than fixed *F* in three
ecosystems out of four but *worse* in eco_r1 (1.34 against 1.01). It protects
the rare species by fishing the abundant ones considerably harder, which moves
the community more overall. The paper reports only the rare-species outcome.

## Phase 2 — the mechanism is not quite the one claimed

Four controls the paper does not run, all on eco1 at the same reference yield:

| rule | worst species vs control | species below 10% | RMS log |
|---|---|---|---|
| fixed *F* ($$\theta = 0$$) | 0.166 | 0 | 0.943 |
| BH<sub>P</sub> | 0.566 | 0 | 0.408 |
| BH<sub>P/B</sub> | 0.061 | 2.6 | 1.352 |
| $$F \propto B$$ | **0.600** | 0 | 0.400 |
| **frozen BH<sub>P</sub> allocation** | **0.411** | 0 | 0.465 |
| frozen BH<sub>P/B</sub> allocation | 0.049 | 3.3 | 1.426 |
| $$\theta = -0.5$$ | ~1e-98 | 5.1 | 15.3 |
| $$\theta = 0.5$$ | 0.579 | 0 | 0.402 |
| $$\theta = 1.5$$ | 0.620 | 0 | 0.412 |

**Production is doing no work.** $$F \propto B_i$$ performs identically to
$$F \propto P_i$$ — marginally better, in fact. Since $$P_i \sim B_i$$ across species, the
active ingredient is proportionality to biomass, and production is a proxy for
it. This matters practically: the paper's own discussion argues for production
partly because it "automatically integrates over all the paths by which biomass
flows into components of an ecosystem", but biomass is very much easier to
estimate, and in this model it does the same job.

**Most of the benefit is the allocation, not the feedback.** Freezing
BH<sub>P</sub>'s year-0 allocation and holding it constant for 50 years reaches
0.411, against 0.566 for the adaptive rule and 0.166 for fixed *F* — about
three-quarters of the benefit on a log scale. The paper's framing is that
fishing "operating adaptively to follow species' production rates over time,
contains a feedback that would help to protect species from overfishing in the
presence of uncertainty". The adaptation does add something real, but simply
*allocating* less fishing mortality to rarer species at the outset, and then
never updating it, recovers most of the gain.

**The benefit saturates, and only the sign of $$\theta$$ is critical.** Writing
the family as $$F_i = c\,B_i^{\theta}$$, the paper's three rules are
$$\theta = 1$$ (BH<sub>P</sub>), $$\theta = 0$$ (fixed) and $$\theta$$ slightly
negative (BH<sub>P/B</sub>). Sweeping $$\theta$$: $$-0.5$$ is catastrophic, $$0$$
is the fixed-*F* baseline, and $$0.5$$, $$1$$ and $$1.5$$ are all much the same.
A harvest control rule needs *some* positive density dependence; how much, past
about $$\theta = 0.5$$, hardly matters.

## Phase 3 item 8b — you cannot switch on larval competition and keep the paper's ecosystem

Appendix B nominates competition among larvae for plankton as the mechanism
that replaces an imposed stock-recruitment relationship. Section 3.6 of the
plan showed that at the published parameter values it is inert. The question
was whether switching it on changes the comparison between rules.

First, a correction to that section. Solving the steady state of Eq. (A.11)
properly, at the published $$r \approx 79\,\mathrm{yr}^{-1}$$ over the larval prey range
the resource level is 1.004 with immigration at its published value and 0.991 with
immigration switched off entirely. **$$I_0 = a_0$$ is not the decisive choice I claimed**; what pins the plankton at carrying capacity is that it regenerates
two orders of magnitude faster than it is grazed. Immigration only matters once $$r_0$$ is already low.

Slowing the plankton with primary production held fixed ($$r_0 a_0$$ constant),
in a single-species screen:

| $$r_0$$ | egg-stage growth elasticity | long-run behaviour |
|---|---|---|
| 10 (published) | −0.009 | stable fixed point |
| 3 | −0.008 | stable fixed point |
| 1 | −0.063 | limit cycle, 1.33–4.63 g m⁻², period ~60 yr |
| 0.5 | −0.152 | limit cycle, 0.46–5.20 g m⁻², period ~60 yr |

There is a Hopf bifurcation between $$r_0 = 3$$ and $$r_0 = 1$$, and larval
competition becomes appreciable at the same place. The cycle period exceeds the
paper's entire 50-year experiment.

Assembling full ecosystems from the same seed as eco1, so that the invader
draws are identical and only the plankton differs:

| | species | $$w_{\max}$$ range | egg elasticity | stable? |
|---|---|---|---|---|
| published ($$r_0 = 10$$) | 15 | 327 g – 13.6 kg | −0.013 | yes |
| $$r_0 = 3$$ | 5 | 112 – 161 g | −0.085 | yes |
| $$r_0 = 1$$ | 6 | 112 – 359 g | −0.125 | yes |

Two things happen at once. Larval competition does switch on, six to ten times
stronger. But **large-bodied species are excluded entirely**: with the plankton
contested, juvenile growth is crushed, and no species reaches even 400 g. The
top third of community biomass sits at 12 g and 1.6 g in the slowed variants,
against 542 g in the published one — so the paper's 400 g fishery would catch
literally nothing.

That is a more fundamental result than the one the experiment was designed to
find. **The paper's ability to assemble a community spanning 100 g to 40 kg,
and therefore to have a large-fish fishery to talk about at all, depends on the
plankton being effectively inexhaustible.** The same assumption that switches
off the recruitment brake is what permits large-bodied life histories.

Two smaller observations. The 60-year limit cycle seen in monoculture does not
survive assembly: all three assembled communities sit at stable fixed points,
so sequential assembly acts as a stability filter, which is what Appendix B
argues it should. And richness falls from 15 to 5–6, with the assembly
algorithm hitting its 40-attempt cap in both slowed variants.

### The rule comparison, once the communities have comparable fisheries

Matching $$w_f$$ on body mass is meaningless when the size structures differ
this much, so it is matched on the fraction of community biomass exposed (48%, the
value the published ecosystem has at 400 g), giving $$w_f$$ of 402 g, 10 g and 1 g. Reference yields then come out comparable (0.198, 0.216, 0.126 g m-2 yr-1),
so the three variants are being fished about equally hard.

| variant | rule | worst species | RMS log |
|---|---|---|---|
| published ($$r_0 = 10$$) | fixed *F* | 0.166 | 0.943 |
| | **BH<sub>P</sub>** | **0.566** | **0.408** |
| | BH<sub>P/B</sub> | 0.061 | 1.352 |
| $$r_0 = 3$$ | fixed *F* | 0.765 | **0.150** |
| | **BH<sub>P</sub>** | **0.940** | 1.029 |
| | BH<sub>P/B</sub> | 0.727 | 0.189 |
| $$r_0 = 1$$ | **fixed *F*** | **0.969** | **0.033** |
| | BH<sub>P</sub> | 0.604 | 1.947 |
| | BH<sub>P/B</sub> | 0.892 | 0.121 |

Two things change, and the second is the headline.

**Fishing stops mattering much at all.** With the plankton contested, the
worst-affected species retains 0.6 to 0.97 of its unfished biomass under every
rule, against 0.06 to 0.57 in the published ecosystem. A real density-dependent
brake protects species from fishing regardless of how the fishing is
distributed, which is exactly what section 3.2 of the plan predicted would
happen if compensatory recruitment were restored by any route.

**The ranking inverts.** At $$r_0 = 1$$, BH<sub>P</sub> leaves its worst species
at 0.604 while fixed *F* leaves it at 0.969, and BH<sub>P</sub> distorts the
assemblage roughly sixty times more on the RMS-log measure (1.95 against 0.03).
The mechanism is not mysterious: $$F \propto B$$ concentrates fishing on the abundant
species, and when every species already has its own brake, that concentration
is pure distortion with nothing to buy. BH<sub>P</sub>'s advantage exists
because the model has no other density dependence for it to be redundant with.

### How much weight this carries

Less than the numbers alone suggest, and the reason is structural rather than
statistical. Contesting the plankton does two things at once: it restores a
recruitment brake, *and* it replaces a 15-species community spanning 327 g to
13.6 kg with a 5- to 6-species community of fish under 360 g. These cannot be
separated within this experiment, because as shown above there is no parameter
setting that gives one without the other. So this is evidence that the
conclusion is contingent on the absence of a recruitment brake, not proof of it.

The clean test is Phase 3 item 8: impose a Beverton-Holt stock-recruitment
relationship on the *published* ecosystem, where the community is held fixed
and only the density dependence changes. That is the experiment to run next, and
this result raises its priority considerably.

## Phase 3 item 8 — the result is contingent on the absence of compensatory recruitment

This is the experiment item 8b could not do. `setBevertonHolt()` imposes a
stock-recruitment relationship on the *published* ecosystem while rescaling
`erepro` so that the initial state is preserved exactly — verified to
$$2\times10^{-4}$$ — so the community, its 15 species, its size structure and
its steady state are **identical across the whole sweep**. Only the sensitivity
of recruitment to egg production changes, from
$$d\ln R_{dd}/d\ln R_{di} = 1$$ at $$L = 0$$ (the paper) to $$0.1$$ at
$$L = 0.9$$.

Worst-affected species, as a fraction of its unfished trajectory, at matched
terminal yield:

| $$L$$ | recruitment elasticity $$1-L$$ | fixed *F* | BH<sub>P</sub> | BH<sub>P/B</sub> | BH<sub>P</sub> advantage |
|---|---|---|---|---|---|
| 0 (the paper) | 1.00 | 0.166 | 0.566 | 0.061 | **3.41×** |
| 0.25 | 0.75 | 0.337 | 0.703 | 0.170 | 2.09× |
| 0.50 | 0.50 | 0.603 | 0.777 | 0.395 | 1.29× |
| 0.75 | 0.25 | 0.778 | 0.824 | 0.635 | 1.06× |
| 0.90 | 0.10 | 0.834 | 0.846 | 0.730 | **1.01×** |

**The advantage closes monotonically and essentially vanishes.** At the paper's
$$L = 0$$, BH<sub>P</sub> leaves its worst species 3.4 times better off than
fixed *F* does. By $$L = 0.5$$ — where halving egg production still costs 29% of
recruitment, which is mild compensation by the standards of fitted
stock-recruitment relationships — the factor is 1.29. By $$L = 0.9$$ it is 1.01.
The same happens to the whole-assemblage measure: the RMS-log ratio between
fixed *F* and BH<sub>P</sub> falls from 2.31 to 1.42. And the species that
BH<sub>P/B</sub> drives below 10% of control at $$L = 0$$ — 2.6 of them — are
all rescued by $$L = 0.25$$ alone, with no change to the harvest rule at all.

So the density dependence that BH<sub>P</sub> supplies through the harvest rule
is largely redundant with the density dependence that a stock-recruitment
relationship supplies through the population. The paper's model has none of the
latter, by deliberate design, and that is why the former looks so valuable.

**This also reconciles item 8b.** There the ranking *inverted*; here it does not
— BH<sub>P</sub> stays nominally best at every level, just by a shrinking
margin. The inversion in item 8b is therefore attributable to the community
change (a 5-species assemblage of fish under 360 g) rather than to the
recruitment brake itself. Two experiments that looked like they disagreed
actually separate cleanly: the brake shrinks the advantage, the community change
reverses it.

**One caveat on the far end of the sweep.** Preserving the steady state requires
$$\mathrm{erepro}(L) = \mathrm{erepro}(0)/(1-L)$$, so $$L = 0.9$$ needs
$$\mathrm{erepro} = 4.0$$, i.e. a reproductive efficiency
$$\epsilon_R = 2.0 > 1$$ — more egg mass than the energy budget allows. That
row is a limiting case, not a plausible parameterisation. The physically
comfortable range is $$L \le 0.5$$ ($$\epsilon_R \le 0.4$$), and the advantage
is already down to 1.29× there.

A side effect worth recording: compensation also damps the model's own
transients. Maximum unfished drift over 50 years falls from 83% at $$L = 0$$ to
21% at $$L = 0.9$$, which is further evidence that the quasi-equilibrium the
paper starts from is as unsettled as it is partly *because* recruitment is
uncompensated.

## Still to do

Phase 3 items 8 (compensatory recruitment), 9 (Eq. A.7), 10–11 ($$w_f$$ and the
P/B measurement range), 12–14 (θ matrix, plankton cap, grid); Phase 4
(update interval and observation error); Phase 5 (replication across 10–20
assemblages). Item 8 remains the most likely single result-changer.
