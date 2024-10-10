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


The "gray area" for extreme values ($inc, \phi \rightarrow 1$ means NaN and/or Inf values are generated when calculating expected testing positivity (prop_pos).  Currently this is due to Qbeta(inverse CDF)=1 for such extreme parameter values. These ranges shrink as `test_prop` increases.

After fixing the pbeta denominator, the "grey area" is away, since we no longer have zero denominator due to overflow.
Now the previosly "grey area" are all in yellow((-10,0) level log differnce) in log-diff_cdf-log_simp.png. This shows clear difference between cdf and log simp method. (simp and log_simp seems merged in this area)

Also, some lwr(qbeta) with output "1" are not -Inf in log_space (e.g. inc=0.3,phi=0.97 or 0.98,t=0.001). but for very extreme cases (e.g. inc=0.3,phi=0.97 or 0.98,t=0.001), the qbeta generate "solid 1" (-Inf in log space), this leads to a numerical jump and makes testing positivity=inc. (should be 1 theoretically, for 0.97,.98 case, it is numerically very close to 1, but suddenly drop due to )
Ongoing: Why? 
Conjecture: qbeta=1 means no one are actually being tested, so the expression of expected positivity fails here

Ongoing: If the log_simp catch the real value better than cdf(pbeta)?? 
