# 2024 Nov 10 (Sun)

We currently have four methods of calculating expected testing positivity for beta distribution based on formula discussed in Machinery of `testing_distrib.rmd`:

-   int: Integrating the beta function ourselves
-   cdf: Using cdf of beta distribution (pbeta function from base R)
-   simp/log: Richard's formula / log version formula (equivalent but log version is slightly better for numerical calculation) as an explicit mathematical simplification
-   est: Hybrid method of log and first order estimation based on cdf formula
    -   For lwr$=Q(1-t)>=1e^-5$ use log method
    -   For lwr$=Q(1-t)<1e^-5$ use first order estimation

All formula is write as function prop_pos_test_new on `testing_funs.R`.

`Expected_Test_positivity_figure.R`: Generate figures for expected test positivity curves as a function of $\phi$ (shape parameter for beta distribution as "test focus"), testing proportion and prevalence.

-   `\pix\test_positivity_vs_phi-test_prop.png`

`Logspace_comparing_methods.R`: Comparing numerical result of methods on log space

-   "cdf" vs "log_simp": `\pix\log-diff_cdf-log_simp.png`

-    "simp" vs" log_simp": `\pix\log-diff_simp-log_simp.png`

Log-space Difference (heat-map color) are function of prevalence (y-axis), $\phi$(x-axis) and testing proportion-test_prop(group). 
For now, the differences are magnitude/absolute value, signs are considered due to log-space calculation, but temporarily ignored. 
Difference of $log(x)-log(y)=0$ cases is highlighted by manually assign some value ($45$ in log-diff_cdf-log_simp.png or $40$ in log-diff_simp-log_simp.png). 
These difference by definition is $-Inf$ but could leads to confusion with the "real gray area" (NaN, caused by numerical flow of methods) in heat-map.

The "gray area" for extreme values (prevalence, $\phi \rightarrow 1$ leads to NaN and/or Inf values when calculating expected testing positivity (prop_pos).
This is due to $Q(1-t)$ which is calculated by Qbeta has numerical overflow and equals to $1$ (inverse CDF)=1 for such extreme parameter values. 
These ranges shrink as `test_prop` increases.

After fixing the denominator, the "grey area" is away, since we no longer have zero denominator due to overflow. 
Now the previosly "grey area" are all in yellow((-10,0) level log differnce) in `\pix\log-diff_cdf-log_simp.png`. 
This shows clear difference between cdf and log simp method. (simp and log_simp seems merged in this area)

Also, some lwr(qbeta) with output "1" are not -Inf in log_space (e.g. inc=0.3,phi=0.97 or 0.98,t=0.001). but for very extreme cases (e.g. inc=0.3,phi=0.97 or 0.98,t=0.001), the qbeta generate "solid 1" (-Inf in log space), this leads to a numerical jump and makes testing positivity=inc. (should be 1 theoretically, for 0.97,.98 case, it is numerically very close to 1, but suddenly drop due to ) Ongoing: Why? Conjecture: qbeta=1 means no one are actually being tested, so the expression of expected positivity fails here

Ongoing: If the log_simp catch the real value better than cdf(pbeta)?? What happend in the yellow area
