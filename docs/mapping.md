# Mapping Law & Plank (2023) onto mizer

Law, R. & Plank, M.J. (2023) *Fishing for biodiversity by balanced harvesting.*
**Fish and Fisheries** 24, 1–17. doi:[10.1111/faf.12705](https://doi.org/10.1111/faf.12705).
Equation numbers below are those of the paper (preprint
doi:[10.1101/2021.06.27.450047](https://doi.org/10.1101/2021.06.27.450047), v3).

Law & Plank write the dynamics in log body mass `x = ln(w/w0)` with `w0 = 1 g`
and state density `u_i(x,t)`; mizer works in `w` with number density
`N_i(w) = phi_i(w)`, and `u = w phi`. Dividing Eq. (A.1) through by `w` turns
it term by term into mizer's equation, so the two are the same PDE.

Units throughout: grams, years, square metres of sea surface. There is no
`kappa`-style rescaling — densities are per m² exactly as in the paper.

## Term by term

| Paper | mizer |
|---|---|
| Eq. (A.1a) somatic growth `-d/dx[eps_i g_i u_i]` | `e_growth`, with `psi = 1 - eps_i` |
| Eq. (A.2) mass-specific intake $$g_i$$ | `alpha * encounter`; `search_vol = gamma w^q` |
| — search-rate coefficient $$A_i$$, exponent $$\alpha$$ | `gamma = A_i/(6 sigma)`, `q = 0.85` |
| — conversion efficiency `K = 0.1` | `alpha = 0.1` (mizer's assimilation efficiency) |
| — no maximum intake rate | `FeedingLevel` overridden to return 0 |
| — no metabolic cost | `ks = k = 0`, so `metab = 0` |
| Eq. (A.3) predation mortality $$d_i$$ | `pred_mort` (standard mizer) |
| Eq. (A.4) box feeding kernel, height `1/(6 sigma)` | `pred_kernel_type = "box"`, `ppmr_min = e^(beta-3 sigma) = 10`, `ppmr_max = e^(beta+3 sigma) = 1e5`; the height is folded into `gamma` because mizer's box kernel has height 1 |
| Eq. (A.6) larval mortality $$\mu_l$$ | `ext_mort` (static in time) |
| Eq. (A.7) background mortality $$\mu_b$$ | added by the custom `Mort` function (food-dependent, see below) |
| Eq. (A.8) allocation $$1-\epsilon_i$$ | `psi = maturity * repro_prop`; `w_mat25 = w_mat/3^(1/rho_m)` gives mizer's maturity exponent `U = rho_m = 15`; `m - n = rho_inf = 0.2` gives `repro_prop = (w/w_max)^0.2` |
| Eq. (A.9) egg production $$R_i$$ | `mizerRDI` with `erepro = 2 eps_R = 0.4` (mizer carries a factor 1/2 for the sex ratio that the paper does not) |
| — no stock-recruitment relationship | `RDD = "noRDD"` |
| Eq. (A.1f) diffusion `(1/2) d/dx(e^-x d/dx[eps_i G_i u_i])` | `use_predation_diffusion = TRUE`; mizer's `D` omits the `eps_i` factor, so the custom `Diffusion` function multiplies it back in |
| Eq. (A.11) plankton logistic + immigration | custom `resource_dynamics = "lpResource"` |
| Eqs (A.12)–(A.14) `r(x)`, `a(x)`, `I(x)` | `resource_rate`, `resource_capacity`, and `other_params$lp_immigration`, each divided by `w` to convert the paper's log-density to mizer's per-mass density |
| Appendix C $$\theta$$ | `interaction` matrix (0.5 on the diagonal, 0.2 off) and `interaction_resource = 1` |
| Eq. (2.3) $$B_i$$ | `lp_biomass()` |
| Eq. (2.4) $$P_i$$ incl. boundary terms | `lp_production()` |
| Eqs (2.7)–(2.9) the three rules | custom `FMort = "lpFMort"` |

## Checks that the mapping is exact

`tests/test_model.R` checks, against hand-coded integrations of the paper's
equations on the model's own grid:

* the encounter rate against Eq. (A.2);
* the predation mortality against Eq. (A.3);
* `psi` against Eq. (A.8);
* the plankton carrying capacity against Eq. (A.13) and the resulting total
  primary production (~4.4e3 g m⁻² yr⁻¹ against the paper's "region of 4000");
* the resource fixed point of `lpResource` against the root of Eq. (A.11);
* `lp_production()` against a finite-difference evaluation of Eq. (2.4).

## The diffusion term

Eq. (A.1f) is written in log mass as

$$
\tfrac12\,\frac{\partial}{\partial x}
  \left( e^{-x}\, \frac{\partial}{\partial x}
  \big[\epsilon_i\, G_i\, u_i\big] \right),
$$

with $$G_i$$ (Eq. A.10) the second moment of prey mass. mizer works in $$w$$ and
adds $$\tfrac12\,\partial^2 (D_i N_i)/\partial w^2$$ with $$D_i$$ in
g² yr⁻¹. The two are the same term. Writing $$B = \epsilon_i G_i u_i$$ and
using $$\partial/\partial x = w\,\partial/\partial w$$ and $$e^{-x} = 1/w$$
(with $$w_0 = 1$$ g),

$$
\tfrac12\,\frac{\partial}{\partial x}
  \left(\frac1w \frac{\partial B}{\partial x}\right)
= \tfrac12\, w\, \frac{\partial}{\partial w}
  \left(\frac1w \cdot w \frac{\partial B}{\partial w}\right)
= \tfrac12\, w\, \frac{\partial^2 B}{\partial w^2},
$$

and the whole of Eq. (A.1) is an equation for $$u = w\phi$$, so dividing
through by $$w$$ leaves $$\tfrac12\,\partial^2 B/\partial w^2$$. Since
$$G_i = D_i/w$$ and $$u_i = w\phi_i$$, we have $$B = \epsilon_i D_i \phi_i$$,
which is mizer's term with the $$\epsilon_i$$ factor that `mizerDiffusion()`
omits and `lpDiffusion()` restores. The same division by $$w$$ turns Eq. (A.1a)
into mizer's advection term, so the correspondence is term-by-term, not just in
aggregate.

## Ambiguities in the paper, and how they are resolved here

1. **Eq. (A.3) as printed evaluates the predator's search area at the prey's
   body mass** — $$A_j (w_0 e^{x})^{\alpha_j}$$, where $$x$$ indexes the prey.
   Eq. (A.2)
   uses the consumer's own size, so this is a typo; mizer uses the predator's
   size and this implementation follows mizer.

2. **Eq. (A.7)** — prose and formula disagree; see [`mu_b.md`](mu_b.md). This
   one changes the results qualitatively, and is exposed as the `mu_b_form`
   option.

3. **Table 1 contradicts the text on $$x_a$$ and $$x_I$$.** The table glosses
   $$w_0 e^{-x_a}$$ as 0.001 g, which with
   $$a(x) = a_0 (w_0 e^{x - x_a})^{1-\lambda}$$ would put the carrying capacity
   at $$2\times10^{9}$$ m⁻² at 1 mg rather than the stated 2000 m⁻². The text's
   $$x_a = \log(0.001)$$ gives $$a(w) = a_0 (w/0.001)^{1-\lambda}$$, which does
   put it at 2000 m⁻² at 1 mg and integrates to a total primary production of
   $$4.4\times10^{3}$$ g m⁻² yr⁻¹ against the paper's "region of 4000". The text
   is therefore taken to be correct.

4. **The invading spectrum.** Appendix B says only "started at a low density on
   a power law, the egg density being 0.002 m^-2". Taken here as the paper's
   own state variable, i.e. the log-density `u(x_0) = 0.002 m^-2`, so
   $$\phi(w_{\mathrm{egg}}) = 2$$ m⁻² g⁻¹, on a power law
   $$u \sim w^{-(\lambda-1)}$$ matching the plankton spectrum. Only the invasion *outcome* depends on this, and only
   weakly: an invader that establishes does so over tens of generations.

## Where this deviates from the paper, and why

1. **The plankton time step.** The paper integrated Eq. (A.11) with explicit
   Euler at `dt = 0.002 yr`. `lpResource()` instead uses the closed-form
   solution of the same Riccati equation over a step, which has the same
   fixed points but is unconditionally stable, so the time step is not tied to
   the fastest-dividing cells (`r` reaches ~320 yr⁻¹ at 1e-10 g).
   `tests/test_convergence.R` checks the two agree to 0.02% at `dt = 0.002`,
   and that `dt = 0.01` and `dt = 0.002` agree to 0.3%. Assembly is run at
   `dt = 0.01`; the finished ecosystems are then relaxed and all fishing runs
   are done at the paper's `dt = 0.002`.

2. **The particular ecosystems differ.** The paper does not report the random
   seed or the realised $$w_{\max}$$ and $$\mu_{\mathrm{egg}}$$ of its 15 species, so the assembly
   procedure of Appendix B is reproduced but not the specific assemblage. The
   figures are therefore qualitatively, not numerically, comparable.

3. **Fig. 2a is empirical data** (an Ecopath model of the West Scotland shelf,
   Alexander et al. 2015) and is not reproduced.
