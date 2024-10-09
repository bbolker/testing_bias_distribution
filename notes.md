e(-10) ## .00004539992976248485

2024 Oct 09 (Wed)
=================

We currently have four methods of calculating expected testing positivity:
* Integrating the beta function ourselves
* Using pbeta
* Richard's formula
* Richard's logspace formula

in prop_pos_test_new.R.
    "cdf"vs"log_simp": log-diff_cdf-log_simp.png
    "simp"vs"log_simp": log-diff_simp-log_simp.png
    Log-space Difference (heatmap color) are function of incidence-inc (y-axis), $\phi$(x-axis) and testing proportion-test_prop(group).
    For now, the differences are magnitude/absoulute value, signs are considered due to logspace calculation, but temporarly ignored.
    Difference of $log(x)-log(y)=0$ cases is highlighted by manually assign some value ($45$ in log-diff_cdf-log_simp.png or $40$ in log-diff_simp-log_simp.png). These difference by definition is -Inf but could leads to confusion with the "real gray area" (NaN) in heatmap.

    The "real gray area" for extreme values ($inc, \phi \rightarrow 1$ means NaN and/or Inf values are genearted when calculating expected testing positivity (prop_pos).  Currently is due to Qbeta(inverse CDF)$=1$ for such extreme parameter values. These range are shrink with test_prop increasing.

