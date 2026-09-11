# The ambiguity in Eq. (A.7), and why this implementation resolves it one way

## The two readings

Appendix A of Law & Plank (2023) introduces background (non-larval,
non-predation) intrinsic mortality like this:

> To set a level playing-field across species, we assume that this is
> proportional to the mass-specific needs for metabolism, **relative to** the
> mass-specific rate at which food becomes available at size `x` [...]. The
> metabolic need should scale with body mass, and we write this as
> `exp(-xi(x - x0,i))` [...]. The mass-specific rate of food intake at size `x`
> relative to size `x0,i` is `g_i(x,t)/g_i(x0,i,t)`.

and then prints

$$
\mu_{b,i}(x,t) \;=\; \mu_{b,i}^{(0)}\,
  \exp\!\big(-\xi\,(x - x_{0,i})\big)\;
  \frac{g_i(x,t)}{g_i(x_{0,i},t)}
  \tag{A.7}
$$

The prose says mortality is metabolic need **divided by** food availability - a
starvation term, large when food is scarce. The printed equation **multiplies**
by the food-availability ratio, so mortality is *small* when food is scarce.
They cannot both be right.

| option | form |
|---|---|
| `mu_b_form = "product"` | $$\mu_b^{(0)}\left(\dfrac{w}{w_0}\right)^{-\xi}\dfrac{g(w,t)}{g(w_0,t)}$$ — as printed |
| `mu_b_form = "ratio"` | $$\mu_b^{(0)}\left(\dfrac{w}{w_0}\right)^{-\xi}\dfrac{g(w_0,t)}{g(w,t)}$$ — as described |

## Why this implementation defaults to "ratio"

The two are not a minor detail: $$g(w,t)/g(w_0,t)$$ falls by about an order of
magnitude from egg size to adult size, so the two forms differ by two orders of
magnitude in adult mortality, and the model behaves completely differently.

Single-species runs to equilibrium, everything else at the paper's values:

| $$w_{\max}$$ | form | mortality at $$w_{\max}/2$$ | total biomass | biomass above 400 g |
|---|---|---|---|---|
| 150 g | product | 0.001 /yr | 10.1 g/m2 | 0 |
| 1 kg | product | 0.000 /yr | 19.8 g/m2 | 19.8 g/m2 |
| 30 kg | product | 0.000 /yr | 34.3 g/m2 | 34.2 g/m2 |
| 150 g | ratio | 0.180 /yr | 3.4 g/m2 | 0 |
| 1 kg | ratio | 0.250 /yr | 3.9 g/m2 | 2.7 g/m2 |
| 30 kg | ratio | 0.601 /yr | 3.9 g/m2 | 3.6 g/m2 |

Three independent anchors in the paper pick out the "ratio" form.

1. **Species biomasses.** Fig. 3a spans about 1e-5 to 1e0 g/m2 across the 15
   species, and Fig. 6 uses the same 1e-6 to 1e0 axis, so the whole assemblage
   totals a few g/m2. Under "product" a *single* species reaches 10-34 g/m2.

2. **Total yield.** The paper reports ~0.25 g m⁻² yr⁻¹ at $$F = 0.1\,\mathrm{yr}^{-1}$$
   over $$w > 400\,\mathrm{g}$$, which needs $$\sum_i B_i(>400\,\mathrm{g}) \approx 2.5$$ g m⁻². Under
   "ratio", one species holds 2.7-3.6 g/m2 above 400 g, so 15 species sharing
   the same plankton production land on that figure. Under "product" the yield
   comes out roughly ten times too large, and the Fogarty ratio would be ~0.6
   per mille rather than the reported 0.06 per mille.

3. **Which species are rare.** In Figs 2c and 3a the small-$$w_{\max}$$ species
   (1, 2, 3, 4) are the common ones and the rare ones sit at larger $$w_{\max}$$.
   Under "product", adult mortality is effectively zero, biomass piles up at
   $$w_{\max}$$ where growth stops, and the largest species always dominates -
   small species cannot invade at all. Under "ratio", mortality rises with $$w_{\max}$$
   (0.18 -> 0.60 /yr in the table above) and the paper's ordering is recovered.

"Product" is therefore taken to be a typo in Eq. (A.7). It remains available:

```r
p <- lp_params(w_max = ..., mu_egg = ..., mu_b_form = "product")
```

Note that "ratio" also makes `mu_b` a genuine starvation feedback - mortality
rises as food per unit metabolic need falls - which is the kind of density
dependence the coexistence argument of Appendix B relies on.
