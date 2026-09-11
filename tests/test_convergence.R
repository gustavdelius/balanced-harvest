# ---------------------------------------------------------------------------
# The two numerical liberties this implementation takes:
#
#   1. the plankton step is the closed-form solution of Eq. (A.11) over a time
#      step rather than the paper's explicit Euler;
#   2. assembly runs at dt = 0.01 yr rather than the paper's dt = 0.002 yr.
#
# Check that neither changes the answer materially.
#
#   Rscript tests/test_convergence.R
# ---------------------------------------------------------------------------

source("R/lp_harvest.R")

failures <- 0
report <- function(label, err, tol) {
    ok <- is.finite(err) && err <= tol
    if (!ok) failures <<- failures + 1
    cat(sprintf("%-4s %-52s max rel diff = %6.2f%%  (tol %.0f%%)\n",
                if (ok) "PASS" else "FAIL", label, 100 * err, 100 * tol))
}

set.seed(11)
w_max  <- exp(runif(6, log(100), log(40000)))
mu_egg <- runif(6, 28, 32)

run <- function(dt, resource_fn = "lpResource", years = 20) {
    p <- lp_params(w_max = w_max, mu_egg = mu_egg)
    resource_dynamics(p) <- resource_fn
    n <- initialN(p)
    for (i in seq_along(w_max)) n[i, ] <- lp_invader_n(p, i)
    initialN(p) <- n
    sim <- project(p, t_max = years, dt = dt, t_save = years,
                   progress_bar = FALSE)
    list(B = getBiomass(sim)[2, ],
         B_pp = sum(NResource(sim)[2, ] * w_full(p) * dw_full(p)))
}

reldiff <- function(a, b) max(abs(a$B / b$B - 1),
                              abs(a$B_pp / b$B_pp - 1))

cat("Running four 20-year projections of a 6-species ecosystem...\n\n")

exact_002 <- run(0.002, "lpResource")
euler_002 <- run(0.002, "lpResourceEuler")
exact_01  <- run(0.01,  "lpResource")

# 1. The closed-form plankton step against the paper's Euler step, both at the
#    paper's dt.  These solve the same ODE, so they must agree.
report("closed-form vs Euler plankton step, dt = 0.002",
       reldiff(exact_002, euler_002), 0.02)

# 2. The assembly time step against the paper's time step.
report("dt = 0.01 vs dt = 0.002", reldiff(exact_01, exact_002), 0.10)

# 3. Euler at dt = 0.01 is outside its stability limit: r reaches ~320/yr at
#    the smallest cell size, so dt*r > 2.  This is why the closed form is used.
r_max <- max(resource_rate(lp_params(w_max[1], mu_egg[1])))
cat(sprintf("\nINFO max plankton intrinsic rate r = %.0f /yr; dt*r = %.2f at dt=0.002, %.2f at dt=0.01\n",
            r_max, 0.002 * r_max, 0.01 * r_max))
cat("     (explicit Euler for the logistic needs dt*r < 2)\n")

cat(if (failures == 0) "\nAll checks passed.\n" else
    sprintf("\n%d CHECK(S) FAILED.\n", failures))
quit(status = if (failures == 0) 0 else 1)
