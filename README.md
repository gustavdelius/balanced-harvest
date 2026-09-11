# Fishing for biodiversity by balanced harvesting, in mizer

A reimplementation of

> Law, R. & Plank, M.J. (2023) *Fishing for biodiversity by balanced
> harvesting.* **Fish and Fisheries** 24, 1–17.
> doi:[10.1111/faf.12705](https://doi.org/10.1111/faf.12705)

using [mizer](https://sizespectrum.org/mizer/).

**Results, figures and comparison with the paper:
<https://gustavdelius.github.io/balanced-harvest/>**

The paper asks whether a moderate amount of fishing, scaled across species in
proportion to each species' production rate, can maintain biodiversity rather
than erode it. It compares three ways of setting the fishing mortality rate of
species *i*:

| | rule | |
|---|---|---|
| fixed | `F_i(t) = F` | Eq. (2.7) |
| BH_P | `F_i(t) = c_P P_i(t)` | Eq. (2.8) |
| BH_P/B | `F_i(t) = c_P/B P_i(t)/B_i(t)` | Eq. (2.9) |

In all three, every fish enters a single mixed-species fishery at `w_f = 400 g`
and everything above that is caught at the same rate, so the balancing is
**across species, not across sizes**. `B_i` and `P_i` are measured over that
same harvested range. `c_P/B` is just a fixed exploitation ratio `E = Y_i/P_i`.

## Layout

```
R/lp_constants.R    every parameter value from Tables 1, 2 and Appendices B, C
R/lp_model.R        the ecosystem model of Eqs (A.1)-(A.14) as a MizerParams
R/lp_assembly.R     sequential assembly of ecosystems (Appendix B)
R/lp_harvest.R      B_i, P_i, the three rules, and yield (Eqs 2.3-2.9)
R/lp_figures.R      Figures 2-6

run_assembly.R      assembles the ecosystems   -> data/ecosystems.rds
run_figures.R       applies the three regimes  -> data/results.rds, figures/

tests/test_model.R        mizer's rates vs hand-coded integrals of the paper
tests/test_convergence.R  the two numerical liberties taken (see docs/mapping.md)
tests/test_ecosystem.R    the assembled ecosystems vs the paper's own criteria

docs/mapping.md     equation-by-equation mapping onto mizer, and the deviations
docs/mu_b.md        an ambiguity in Eq. (A.7) that changes the results
docs/robustness.md  which of the paper's arbitrary choices could matter, and a
                    plan for testing them (design document, not yet run)
```

## Running it

```bash
Rscript tests/test_model.R        # no inputs needed
Rscript run_assembly.R            # -> data/ecosystems.rds
Rscript tests/test_ecosystem.R    # needs data/ecosystems.rds
Rscript run_figures.R             # -> data/results.rds, figures/
Rscript tests/test_convergence.R  # slow; no inputs needed
```

`run_assembly.R` takes tens of minutes: it assembles four ecosystems, each by
up to 40 sequential invasions with a 50-year relaxation after every one.

## What the model is

Everything is in grams, years and square metres of sea surface — the paper's
own units, with no `kappa`-style rescaling. Fifteen fish species share a life
history template and differ only in maximum body mass (log-uniform over
100 g – 40 kg) and larval death rate, and are assembled one at a time onto a
plankton community spanning 1e-10 g to 1 g.

Points where it differs from a stock mizer model:

* **no satiation and no metabolic cost** — energy available is `K` times the
  encounter rate, so the feeding level is identically zero;
* **no stock-recruitment relationship** — `RDD = noRDD`, all density dependence
  is emergent;
* **intrinsic mortality is food-dependent** — a larval term plus a background
  term that responds to the current growth rate, so it needs a custom `Mort`;
* **the plankton follow a logistic with immigration**, not mizer's
  semichemostat, so they need custom resource dynamics;
* **the fishing rules are a feedback on the model's own state**, so they need a
  custom `FMort` — `project()`'s `effort` argument has to be prescribed in
  advance and cannot see the evolving ecosystem.

`docs/mapping.md` gives the full equation-by-equation correspondence, and
`tests/test_model.R` checks it numerically against the paper's equations.

## Two things to know before trusting the numbers

1. **Eq. (A.7) is ambiguous** — its prose and its printed formula disagree, and
   the choice changes the results qualitatively. See [`docs/mu_b.md`](docs/mu_b.md)
   for the evidence used to resolve it.

2. **The ecosystems are not the paper's ecosystems.** Law & Plank do not report
   the random seed or the realised life histories of their 15 species, so this
   reproduces the assembly *procedure* of Appendix B, not the particular
   assemblage. The figures are comparable in kind, not number for number.

## Unrelated files

`bh_model.R`, `bh_run.R`, `bh_figures.R`, `bh_rare*.R` and `bh_results.rds` are
an earlier, separate exercise that applied the same three fishing rules to
mizer's built-in North Sea model (`NS_params`). They are not part of this
reimplementation and do not use the paper's ecosystem model.
