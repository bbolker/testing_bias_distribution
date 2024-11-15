---
output: html_document
editor_options: 
  chunk_output_type: inline
---

# 2024 Nov 10 (Sun)

## Beta Distribution

`testing_distrib.rmd` general explanation of motivation of the project. With methodology of beta distribution idea and some preliminary simluation.

`Expected_Test_positivity_figure.R`: Generate figures for expected test positivity curves as a function of $\phi$ (shape parameter for beta distribution as "test focus"), testing proportion and prevalence.

-   `\pix\test_positivity_vs_phi-test_prop.png`: group by different test proportion

-   `\pix\test_positivity_vs_test_proportion-phi.png`: group by different $\phi$ value

### Numerical issues

#### qbeta overflow

We currently have four methods of calculating expected testing positivity for beta distribution based on formula discussed in Machinery of `testing_distrib.rmd`:

-   int: Integrating the beta function ourselves $$
    \frac{\left(\int_{Q(1-t)}^1 B(y,\phi) y \, dy\right)}{1-\Phi_B(Q(1-t))} = 
    \frac{\left(\int_{Q(1-t)}^1 B(y,\phi) y \, dy\right)}{t}
    $$

-   cdf: Using cdf of beta distribution (pbeta function from base R) $$
    \frac{CDF_{beta}(Q(1-t),a+1,b)}{t}
    $$ where $a,b$ are shape parameter after parameterization by prevalence $i$ and $\phi$. $t$ is test proportion.

-   simp/log: Richard's formula / log version formula (equivalent but log version is slightly better for numerical calculation) as an explicit mathematical simplification of cdf method $$
    i(1+\frac{Q(1-t)^a(1-Q(1-t))^b}{ta\mathbb{B}(a,b)})
    $$\

-   est: Hybrid method of log and first order estimation based on cdf formula

    -   For lwr$=Q(1-t)>=1e^{-5}$ use log method
    -   For lwr$=Q(1-t)<1e^{-5}$ use first order estimation $$
        \frac{b+1-(1-Q(1-t))ab}{b+1-(1-Q(1-t))(1-a)b}
        $$

All formula is write as function prop_pos_test_new on `testing_funs.R`.

`Logspace_comparing_methods.R`: Comparing numerical result of methods on log space. Log-space Difference (heat-map color) are function of prevalence (y-axis), $\phi$(x-axis) and testing proportion-test_prop(group). For now, the differences are magnitude/absolute value, signs are considered due to log-space calculation, but temporarily ignored.

Difference of $log(x)-log(y)=0$ cases is highlighted by manually assign some value ($45$ in `\pix\log-diff_cdf-log_simp.png` or $40$ in `\pix\log-diff_simp-log_simp.png`). These difference by definition is $-Inf$ but could leads to confusion with the "real gray area" (NaN, caused by numerical flow of methods) in heat-map.

The "gray area" for extreme values (prevalence, $\phi \rightarrow 1$ leads to NaN and/or Inf values when calculating expected testing positivity (prop_pos). This is due to $Q(1-t)$ which is calculated by Qbeta has numerical overflow and equals to $1$ (inverse CDF)=1 for such extreme parameter values. These ranges shrink as `test_prop` increases.

-   "cdf" vs "log_simp": `\pix\log-diff_cdf-log_simp.png`

After fixing the denominator issue (Qbeta appears in denominator, now is replaced by testing proportion), the "grey area" is away, since we no longer have zero denominator due to overflow. Now the previously "grey area" are all in yellow((-10,0) level log difference) in `\pix\log-diff_cdf-log_simp.png`. This shows clear difference between cdf and log simp method. (simp and log_simp seems merged in this area)

-   "simp" vs "log_simp": `\pix\log-diff_simp-log_simp.png`

log_simp method works well when away from the "grey area". Numerical overflow $Q(1-t)$ is still causing issues near and in "grey area".

-   (TO DO?) one way to "resolve" this is by increase numerical accuracy of qbeta function in R. qbeta use a newton-like methods and limit the maximum iteration to 1000. Idea: we could try to recreate the algorithm in r level (available on Github) and then maybe try C level.

-   "est": use log_simp when $Q(1-t)$ is away from $1$ (away from "grey area") and use first order estimation of "cdf" formula when it is close to 1

    -   Current boundary is set to $1e^{-5}$ (TO DO? better analysis). where the the estimation and log_simp are close enough, but before the log_simp starting to break down
    -   "est" vs "log_simp" for whole parameter space: `\pix\log-diff_est-log.png`
    -   `\pix\OLD-log-diff_est-log.png` & `\pix\New-log-diff_est-log.png`: A focuse on a certain test proportion. OLD is estimation directly compare with log_simp. New is the hybrid "est" method compare with logsimp, grey area is where $Q>=1^e{-5}$ both method use log_simp formula. This is to show the boundary is reasonable
    -   `\pix\slice_log_diff_est-logsimp.png`: A more specific slice at $\phi=0.851$.

#### Qbeta vs qbeta

Ben previously have some other issue for qbeta some times generate error, so a Qbeta function is created in `testing_funs.R` to avoid error.

`Qbeta_Issues.R` is trying to see if qbeta will causing us such problems, thus if Qbeta is necessary.

-   `\pix\log-diff_Qbeta-qbeta.png`: seems not.
