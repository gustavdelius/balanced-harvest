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
R/lp_constants.R          every parameter value from Tables 1, 2 and
                          Appendices B, C
R/lp_model.R              the ecosystem model of Eqs (A.1)-(A.14) as a
                          MizerParams object
R/lp_assembly.R           sequential assembly of ecosystems (Appendix B)
R/lp_harvest.R            B_i, P_i, the three rules, and yield (Eqs 2.3-2.9)
R/lp_figures.R            Figures 2-6
R/lp_experiments.R        the yield-biodiversity frontier machinery the
                          robustness investigation runs on
R/lp_resource_variants.R  larval competition for plankton, switched on
R/ns_*.R                  the separate North Sea exercise: survey biomass,
                          FishBase species parameters, interaction matrix,
                          model, runs and figures

run_assembly.R            the four ecosystems     -> data/ecosystems.rds
run_figures.R             the three regimes       -> data/results.rds,
                                                     figures/
run_robustness.R          phases 1 and 2          -> data/robustness/
run_recruitment.R         phase 3 item 8          -> data/recruitment/
run_resource_sweep.R      phase 3 item 8b         -> data/resource/
run_fishery_design.R      phase 3 items 10, 11    -> data/fishery/
run_implementation.R      phase 4                 -> data/implementation/
run_replication.R         phase 5, 13 ecosystems  -> data/replication/
run_egg_assembly.R        egg mass as an axis     -> data/egg/
run_activity_assembly.R   a fast-slow continuum   -> data/activity/
run_activity_harvest.R    the rules on it         -> data/activity_harvest/
run_activity_bhb.R        BH_B added to those     -> data/activity_harvest/
R/ns_pdf.R                the North Sea figures   -> balanced_harvesting.pdf

tests/test_model.R        mizer's rates vs hand-coded integrals of the paper
tests/test_convergence.R  the two numerical liberties taken (see
                          docs/mapping.md)
tests/test_ecosystem.R    the assembled ecosystems vs the paper's own criteria

docs/index.md             the paper reproduced in mizer, figure by figure
docs/mapping.md           equation-by-equation mapping onto mizer, and the
                          deviations
docs/mu_b.md              an ambiguity in Eq. (A.7) that changes the results
docs/robustness.md        which of the paper's arbitrary choices could matter,
                          and a plan for testing them
docs/robustness-results.md
                          what the experiments found
docs/life-history.md      where B ~ P comes from, and two attempts to break it
docs/north-sea.md         the three rules on mizer's North Sea model
docs/pdf/                 a PDF of each of those pages

make_pdfs.py              renders every page listed in docs/_data/nav.yml to
                          a PDF in docs/pdf/
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

That is the paper itself. The investigation on top of it is driven by the other
`run_*.R` scripts, each caching into its own directory under `data/` so that a
rerun picks up where it left off; `docs/robustness.md` says what each phase
asks and `docs/robustness-results.md` what it found. The four that build the
alternative ecosystems take arguments:

```bash
Rscript run_egg_assembly.R      <seed> <varied|fixed>
Rscript run_activity_assembly.R <seed> <neutral|gradient>
Rscript run_activity_harvest.R  <seed>
Rscript run_activity_bhb.R      <seed>
```

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

## A second test, on a real ecosystem

[`docs/north-sea.md`](docs/north-sea.md) applies the same three rules to mizer's
built-in North Sea model, with real depleted North Sea species (thornback ray,
Atlantic halibut) added at abundances derived from ICES IBTS survey data. It is
a separate exercise — it does not use the paper's ecosystem model — and it asks
a different question: does the conclusion survive contact with a calibrated
ecosystem and measured abundances?

Short answer: yes, but only once the ecosystem contains something capable of
being lost, and the usual "most depleted species" diagnostic points the wrong
way. Code in `R/ns_*.R`.

## The site as PDFs

Every page of <https://gustavdelius.github.io/balanced-harvest/> is also kept as
a PDF in `docs/pdf/`, reachable from the `PDF` link in the site's header bar.
Regenerate them after editing any page:

```bash
python3 make_pdfs.py
```

It needs `pandoc`, `google-chrome`, `pdflatex` and Python's `pypdf`, but not
Jekyll: pandoc turns each page into print-styled HTML with the maths as MathML,
headless Chrome prints it, and a small LaTeX overlay adds the page numbers.

The list of pages comes from `docs/_data/nav.yml`, so a page added to the site
navigation gets a PDF without the script being touched, and the header link in
`docs/_layouts/default.html` finds it by name (`/` to `index.pdf`,
`/name.html` to `name.pdf`). The one convention the pages have to keep is in
the maths: it is written the way kramdown needs it, `$$...$$` for inline as
well as display, so the script takes a `$$` alone on its own line as a display
delimiter and treats every other `$$` as inline. Do not start a display block
with the formula on the same line as its `$$`.
