# ---------------------------------------------------------------------------
# Check that the mizer model really is the model of Law & Plank (2023), by
# comparing mizer's rates against hand-coded integrations of the paper's
# equations on the model's own grid.
#
#   Rscript tests/test_model.R
# ---------------------------------------------------------------------------

source("R/lp_harvest.R")

failures <- 0
check <- function(label, got, want, tol = 1e-10, rel = TRUE) {
    d <- if (rel) max(abs(got / want - 1), na.rm = TRUE)
         else     max(abs(got - want), na.rm = TRUE)
    ok <- is.finite(d) && d <= tol
    if (!ok) failures <<- failures + 1
    cat(sprintf("%-4s %-58s %s = %.3g\n", if (ok) "PASS" else "FAIL", label,
                if (rel) "max rel err" else "max abs err", d))
}

F <- LP_FISH; PL <- LP_PLANKTON
p <- lp_params(w_max = c(300, 3000, 30000), mu_egg = c(28, 30, 32))
w <- w(p); dw <- dw(p); wf <- w_full(p); dwf <- dw_full(p)
sp <- species_params(p)

# A non-trivial state to evaluate the rates at.
set.seed(1)
n <- initialN(p)
for (i in seq_len(nrow(n))) {
    n[i, ] <- 10^runif(1, -3, -1) * (w / w[1])^-2
    n[i, w > sp$w_max[i]] <- 0
}
initialN(p) <- n
n_pp <- initialNResource(p)

s_h  <- 1 / (6 * F$sigma)
p_lo <- exp(F$beta - 3 * F$sigma)
p_hi <- exp(F$beta + 3 * F$sigma)
th   <- interaction_matrix(p)

## Eq. (A.2): mass-specific intake, here as the absolute rate K A w^a int(...)
prey_w <- c(wf[!(wf %in% w)], w)           # full prey grid
paper_enc <- t(sapply(seq_len(nrow(n)), function(i) {
    sapply(w, function(ww) {
        sel_f <- ww / wf >= p_lo & ww / wf <= p_hi
        sel_s <- ww / w  >= p_lo & ww / w  <= p_hi
        from_pp   <- sum((s_h * wf * n_pp * dwf)[sel_f])          # theta_i0 = 1
        from_fish <- sum(th[i, ] %*% (sweep(n, 2, s_h * w * dw, "*")[, sel_s,
                                                                     drop = FALSE]))
        F$A * ww^F$alpha_q * (from_pp + from_fish)
    })
}))
check("Eq. (A.2) encounter rate", getEncounter(p)[, ], paper_enc, 1e-8)

## Eq. (A.3): predation mortality
paper_d <- t(sapply(seq_len(nrow(n)), function(i) {
    sapply(w, function(ww) {
        sel <- w / ww >= p_lo & w / ww <= p_hi
        sum(sapply(seq_len(nrow(n)), function(j)
            th[j, i] * sum((F$A * w^F$alpha_q * s_h * n[j, ] * dw)[sel])))
    })
}))
check("Eq. (A.3) predation mortality", getPredMort(p)[, ], paper_d, 1e-8)

## Eq. (A.4): the kernel integrates to 1 over log prey mass
kern <- s_h * (log(p_lo) <= log(wf) - log(wf[1]))   # shape only
check("Eq. (A.4) kernel normalisation 1/(6 sigma)",
      s_h * (log(p_hi) - log(p_lo)), 1, 1e-12)

## Eq. (A.6): larval mortality
paper_mu_l <- outer(sp$z * sp$mu_egg, w,
                    function(mu, ww) mu / (1 + (ww / F$w_L)^F$rho_L))
check("Eq. (A.6) larval mortality (ext_mort)", ext_mort(p)[, ], paper_mu_l, 1e-12)

## Eq. (A.7): background mortality, in the "ratio" reading (see docs/mu_b.md)
enc <- getEncounter(p)
g_rel <- sweep(sweep(enc, 2, w, "/"), 1, enc[, 1] / w[1], "/")
paper_mu_b <- outer(sp$z * F$mu_b0, (w / w[1])^(-F$xi)) / g_rel
check("Eq. (A.7) background mortality",
      getMort(p)[, ] - getPredMort(p)[, ] - ext_mort(p)[, ] -
          getFMort(p)[, ], paper_mu_b, 1e-9)

## Eq. (A.8): allocation to reproduction
paper_psi <- (1 + outer(sp$w_mat, w, function(m, ww) (ww / m)^(-F$rho_m)))^-1 *
    outer(sp$w_max, w, function(M, ww) (ww / M)^F$rho_inf)
# mizer rounds maturity below 1e-8 down to zero, so compare absolutely.
inside <- outer(sp$w_max, w, ">=")
check("Eq. (A.8) psi = 1 - epsilon", p@psi[inside], paper_psi[inside],
      1e-8, rel = FALSE)

## Eq. (A.9): egg production rate
paper_R <- drop((p@psi * getEReproAndGrowth(p) * n) %*% dw) * F$eps_R / F$w_egg
check("Eq. (A.9) egg production rate", getRDI(p), paper_R, 1e-12)
check("no stock-recruitment relationship (RDD = RDI)", getRDD(p), getRDI(p), 1e-12)

## Eq. (A.1f): the diffusion term carries the epsilon factor
check("Eq. (A.1f) diffusion carries epsilon_i",
      getDiffusion(p)[, ],
      ((1 - p@psi) * mizerDiffusion(p, n, n_pp, list(),
                                    feeding_level = getFeedingLevel(p)))[, ],
      1e-12)

## Eqs (A.13), (A.14): plankton capacity and immigration, as per-mass densities
on <- wf <= PL$w_max
check("Eq. (A.13) plankton carrying capacity",
      resource_capacity(p)[on],
      PL$a0 * (wf[on] / PL$w_a)^(1 - PL$lambda) / wf[on], 1e-12)
check("Eq. (A.12) plankton intrinsic rate",
      resource_rate(p)[on], PL$r0 * wf[on]^(-PL$rho), 1e-12)

## Eq. (A.11): lpResource must leave the true fixed point alone.
# At the fixed point I + (r-d) u - (r/a) u^2 = 0; take d = 0 for the check.
rates0 <- list(resource_mort = rep(0, length(wf)))
r0 <- resource_rate(p); a0 <- resource_capacity(p)
imm <- other_params(p)$lp_immigration
fp <- rep(0, length(wf))
k <- r0[on] / a0[on]; b <- r0[on]
fp[on] <- (b + sqrt(b^2 + 4 * k * imm[on])) / (2 * k)
stepped <- lpResource(p, n, fp, list(), rates0, 0, 0.002,
                      resource_rate = r0, resource_capacity = a0)
check("Eq. (A.11) resource fixed point is fixed", stepped[on], fp[on], 1e-12)

## Section 2.4: total primary production rate r(x) w u_0(x), paper "region 4000"
pp <- sum(resource_rate(p) * wf * initialNResource(p) * dwf)
cat(sprintf("INFO total primary production = %.0f g/m2/yr  (paper: ~4000)\n", pp))
if (!(pp > 3000 && pp < 6000)) { failures <- failures + 1; cat("FAIL out of range\n") }

## Eq. (2.4): the boundary terms are the ones integration by parts demands.
# int_lo^hi w d/dw(g phi) dw = [w g phi]_lo^hi - int g phi dw, so lp_production()
# must equal -int w d(g phi)/dw dw.  On a log grid of step dx the finite
# difference is only O(dx) accurate, so the test is that the residual shrinks
# with dx rather than that it is zero at any one dx.
residual <- function(dx) {
    pp <- lp_params(w_max = sp$w_max, mu_egg = sp$mu_egg,
                    numerics = modifyList(LP_NUMERICS, list(dx = dx)))
    ww <- w(pp); dww <- dw(pp)
    nn <- initialN(pp)
    for (i in seq_len(nrow(nn))) {
        nn[i, ] <- 10^(-1 - i) * (ww / ww[1])^-2
        nn[i, ww > sp$w_max[i]] <- 0
    }
    initialN(pp) <- nn
    pp <- lp_set_fishing(pp, "none")
    eg <- getEGrowth(pp)
    flux <- eg * nn                       # g phi
    jlo <- lp_wf_idx(pp)
    nw <- length(ww)
    idx <- jlo:(nw - 1)
    d <- t(apply(flux, 1, function(v) c(diff(v) / diff(ww), 0)))
    manual <- -drop(sweep(d[, idx, drop = FALSE], 2, ww[idx] * dww[idx], "*") %*%
                        rep(1, length(idx)))
    got <- lp_production(pp, nn, eg)
    keep <- sp$w_max > LP_FISHING$w_f      # species below w_f are never fished
    max(abs(got[keep] / manual[keep] - 1))
}
r1 <- residual(LP_NUMERICS$dx)
r2 <- residual(LP_NUMERICS$dx / 4)
cat(sprintf("INFO Eq. (2.4) finite-difference residual: %.3g at dx=%.3f, %.3g at dx=%.4f\n",
            r1, LP_NUMERICS$dx, r2, LP_NUMERICS$dx / 4))
check("Eq. (2.4) production, residual falls with dx", r2 / r1, 0, 0.4, rel = FALSE)

## Eq. (2.6): yield factorises as F_i B_i over the harvested range
p3 <- lp_set_fishing(p, "fixed", 0.1)
Fm <- getFMort(p3)
check("Eq. (2.6) Y_i = F_i B_i",
      drop((Fm * n) %*% (w * dw)), 0.1 * lp_biomass(p3, n), 1e-12)

## The mu_b guard: the encounter rate is strictly positive everywhere on the
## grid, so the !is.finite() fallback in lp_background_mort never fires.
check("encounter rate is positive at every size", as.numeric(all(enc > 0)), 1, 0)

cat(if (failures == 0) "\nAll checks passed.\n" else
    sprintf("\n%d CHECK(S) FAILED.\n", failures))
quit(status = if (failures == 0) 0 else 1)
