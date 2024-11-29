---
title: "Notes"
author: "Richard Zhao"
date: "2024-11-26"
knit: (function(inputFile, encoding) {
  rmarkdown::render(inputFile, encoding = encoding, output_dir = "outputs") })
output:
  html_document:
    keep_md: true
---



# 2024 Nov 10 (Sun)

Notes grouped by project:

## Beta Distribution

[`testing_distrib.rmd`](testing_distrib.rmd) is the general explanation of motivation of the project. With methodology of beta distribution idea and some preliminary simulation.

[`Expected_Test_positivity_figure.R`](Expected_Test_positivity_figure.R): Generate figures for expected test positivity curves as a function of $\phi$ (shape parameter for beta distribution as "test focus"), testing proportion and prevalence.

-   [`./pix/test_positivity_vs_phi-test_prop.png`](./pix/test_positivity_vs_phi-test_prop.png): group by different test proportion

-   [`./pix/test_positivity_vs_test_proportion-phi.png`](./pix/test_positivity_vs_test_proportion-phi.png): group by different $\phi$ value

    -   @bbolker 's new color scaling applied

### Ratio of prevalence and test positivity(pos_prop)

Part of motivation for the project is find relationship into relationship between prevalence(prev) and test positivity(pos_prop)

[`inc-testing_positivity-ratio.R`](inc-testing_positivity-ratio.R): Try to dig further with the ratio of prev/pos_prop under different $\phi$ and test proportion(test_prop).

-   [`./pix/prev-pos_prop-ratio_prev_slice.png`](./pix/prev-pos_prop-ratio_prev_slice.png):prev/pos_prop ratio as a function of $\phi$, prev and test_prop, x-axis is $\phi$, grouped by different prev, and colored by test_proportion.

-   [`./pix/prev-pos_prop-ratio_test_prop_slice.png`](./pix/prev-pos_prop-ratio_test_prop_slice.png):prev/pos_prop ratio as a function of $\phi$, prev and test_prop, x-axis is test proportion, grouped by different prev, and colored by $\phi$.

    -   (TO DO) what happened at the jump/cusp?? Another over flow issue for $Q(1-t)$ is too close to $0$!!!!

### Methods of calculate expected testing positivity

We currently have four methods of calculating expected testing positivity(pos_prop) $\bar{p}$ for beta distribution based on formula discussed in Machinery of `testing_distrib.rmd`:

-   int: Integrating the beta function ourselves $$
    \bar{p}=\frac{\left(\int_{Q(1-t)}^1 B(y,\phi) y \, dy\right)}{1-\Phi_B(Q(1-t))} = 
    \frac{\left(\int_{Q(1-t)}^1 B(y,\phi) y \, dy\right)}{t}
    $$

-   cdf: Using cdf of beta distribution (pbeta function from base R) $$
    \bar{p}=\frac{CDF_{\beta}(Q(1-t),a+1,b)}{t}
    $$ where $a,b$ are shape parameter after parameterization by prevalence $i$ and $\phi$. $t$ is test proportion.

-   simp/log: Richard's formula / log version formula (equivalent but log version is slightly better for numerical calculation) as an explicit mathematical simplification of cdf method $$
    \bar{p}=i(1+\frac{Q(1-t)^a(1-Q(1-t))^b}{t a\mathbb{B}(a,b)})
    $$ where $\mathbb{B}(a,b)$ is the complete beta function.

-   est: Hybrid method of log and first order estimation based on cdf formula

    -   For $1-lwr=1-Q(1-t)>=1e^{-5}$ use log method
    -   For $1-lwr=1-Q(1-t)<1e^{-5}$ use first order estimation $$
        \bar{p} \approx \frac{b+1-(1-Q(1-t))ab}{b+1-(1-Q(1-t))(1-a)b}
        $$

All formula is write as function prop_pos_test_new on `testing_funs.R`.

[`Logspace_comparing_methods.R`](Logspace_comparing_methods.R): Comparing numerical result of methods on log space. Log-space Difference (heat-map color) are function of prevalence (y-axis), $\phi$(x-axis) and testing proportion-test_prop(group). For now, the differences are magnitude/absolute value, signs are considered due to log-space calculation, but temporarily ignored.

Difference of $log(x)-log(y)=0$ cases is highlighted by manually assign some value ($45$ in [`./pix/log-diff_cdf-log_simp.png`](./pix/log-diff_cdf-log_simp.png) or $40$ in [`./pix/log-diff_simp-log_simp.png`](./pix/log-diff_simp-log_simp.png)). These difference by definition is -Inf but could leads to confusion with the "real gray area" (NaN, caused by numerical flow of methods) in heat-map.

[`betaParams.md`](betaParams.md): (TO DO?) @dushoff 's suggestion to consider other parameterization of beta distribution shape parameter.

### Numerical issues

#### qbeta overflow

The "gray area" for extreme values (prevalence, $\phi \rightarrow 1$ leads to NaN and/or Inf values when calculating expected testing positivity (prop_pos). This is due to $Q(1-t)$ which is calculated by qbeta has numerical overflow and incorrectly gives $1$ for such extreme parameter values. These ranges shrink as test_prop increases.

-   [`./pix/log-diff_cdf-log_simp.png`](./pix/log-diff_cdf-log_simp.png) :"cdf" vs "log_simp":
-   [`./pix/log-diff_simp-log.png`](./pix/log-diff_simp-log.png) :"simp" vs "log_simp"

After fixing the denominator issue (qbeta appears in denominator, now is replaced by testing proportion), the "grey area" is away, since we no longer have zero denominator due to overflow. This shows difference between cdf and log simp method. (simp and log_simp seems merged). Now the previously "grey area" are not NaN values and all in yellow((-10,0) level log difference) in [`./pix/log-diff_simp-log.png`](./pix/log-diff_simp-log.png).

Still, in this "grey/yellow area", where numerical overflow leads to $Q(1-t)=1$ incorrectly, logsimp method gives $\bar{p}=i$ incorrectly. Near "grey/yellow area", log_simp might also generate numerical value $>1$ for $\bar{p}$. log_simp method works well when away from the "grey area".

-   (TO DO?) one way to "resolve" this is increase numerical accuracy of qbeta function in R. qbeta use a newton-like methods and limit the maximum iteration to 1000. Idea: we could try to recreate the algorithm in r level (available on Github) and then maybe try C level.

-   Another way "est": use log_simp when $Q(1-t)$ is away from $1$ (away from "grey area") and use first order estimation of "cdf" formula when it is close to 1

    -   Current boundary is set to $1-Q(1-t)=1e^{-5}$ (TO DO? better analysis). where the the estimation and log_simp are close enough, but before the log_simp starting to break down
    -   ???Also a boundary for $Q(1-t)=1e^{-5}$ to avoid similar problem of $Q(1-t)=0$ makes overflow too

-   [`./pix/log-diff_est-log.png`](./pix/log-diff_est-log.png): "est" vs "log_simp" for whole parameter space:

-   [`./pix/OLD-log-diff_est-log.png`](./pix/OLD-log-diff_est-log.png): Comparison directly between log_simp and first order estimation, focused on a certain test proportion.

-   [`./pix/New-log-diff_est-log.png`](./pix/New-log-diff_est-log.png): . This is comparison between hybrid "est" method and log_simp, grey area (log(0)) is where $Q>=1^e{-5}$ and both method use log_simp formula. This is to show the boundary $1e^{-5}$ is reasonable

-   [`./pix/slice_log_diff_est-logsimp.png`](./pix/slice_log_diff_est-logsimp.png): A more specific slice at $\phi=0.851$. The overflow problem of log_simp method is more obvious.

#### Qbeta vs qbeta

@bbolker previously have some other issue for qbeta some times generate error, so a Qbeta function is created in `testing_funs.R` to avoid error.

[`Qbeta_Issues.R`](Qbeta_Issues.R) is trying to see if qbeta will causing us such problems, thus if Qbeta is necessary.

-   [`./pix/log-diff_Qbeta-qbeta.png`](./pix/log-diff_Qbeta-qbeta.png): qbeta and Qbeta agree with each other on every point
-   [`./pix/qbeta_logsimp.png`](./pix/qbeta_logsimp.png): qbeta works fine with log_simp method

## Hazard Ratio

### Assumptions and notations

Assuming at any time $t$,

-   Every individual in population have the same baseline cumulative hazard $b$ of being tested, which is irrelevant to being infected or not.

    -   (??) if $b$ corresponding to the behavior factors should be a variable $b(t)$ of $t$?

-   Assuming infection provides a (constant) hazard offset $\phi$ for being tested.

-   $T$ is the testing proportion of the whole population (also the average testing probability).

-   $P$ is the testing positivity.

-   $Y$ is the prevalence (also the average infection probability)

### Target

We would like to use time series of $T$ and $P$ as observed from data and trying to inference $Y$ and $\phi$ with a statistical framework.

### Model

Then the probability/risk of an uninfected individual being tested is governed by the baseline hazard $b$: $$ T_B=r(b)=1-e^{-b} $$And the probability/risk of an infected individual being tested is governed by the baseline hazard $b$ and the offset $\phi$: $$ T_Y=r(b+\phi)=1-e^{-(b+\phi)} $$

For simplicity we denote $B=e^{-b}$ and $\Phi=e^{-\phi}$, then we have $$ T_B=1-B $$$$ T_Y=1-B\Phi$$Then at any time, the testing proportion $T$, also as the probability of a random individual being tested, can be expressed as: $$ T=T_B (1-Y)+T_Y Y=(1-Y)(1-B)+Y(1-B\Phi)$$This gives us $B$ (thus $b=-log(B)$ as the baseline hazard) $$B(Y,T,\Phi)=\frac{1-T}{1-Y(1-\Phi)}$$

-   Constraint: $B\leq 1$ for combination of $T,Y,\Phi$:
    -   For $B\leq 1$, $\Phi$ have a constraint $\Phi \geq 1-T/Y$.
    -   If $T \geq Y$, all feasible value of $\Phi \in [0,1]$ works.
    -   If $T < Y$ (which is more possible in reality), we must have $\Phi \geq 1-T/Y$ as a lower bound $\Leftrightarrow \phi$ has a upper bound.

Consider we can observe testing positivity $P$ and testing proportion $T$ from data. Try to interpret previous function for $B$ (thus $b$) as a function depend on $T$, $Y$, $\Phi$.

**(To be discussed)** $B$ can be seen as the behavioral(???) factor reflect public and administrative reaction to the endemic. They partially observe $Y$, and the testing probability is also limited by testing availability reflected by $T$. While $\Phi$ (thus $\phi$) is more about medical/biological factor of the infection and testing strategy, like the severity/clarity of the symptoms and if the test focus on high risk people

**(To be discussed)** Interpretation of $B(T,Y,\Phi)$:

-   $T \uparrow \Rightarrow B \downarrow \Rightarrow T_B=1-B \uparrow \ \& \ T_Y=1-B\Phi \uparrow$: more testing available, more people get tested in general.

-   (???) $Y \uparrow \Rightarrow B \uparrow \Rightarrow T_B=1-B \downarrow \ \& \ T_Y=1-B\Phi \downarrow$: what does this catch?

-   $\phi \uparrow \Rightarrow \Phi=e^{-\phi} \downarrow \Rightarrow B \uparrow$ due to decrease denominator:

    -   $T_B=1-B \downarrow$: the symptom of the disease is more clear/severe or the test is targeting more on high risk population, thus uninfected people is less likely of being tested

    -   $T_Y=1-B\Phi ??=1-\frac{(1-T)\Phi}{1-Y(1-\Phi)}$: infected people are more likely to get tested (Not necessary???)

Current simplify idea is assume $\phi$ is a constant for a certain period of time and in future it could be a function of $t$.

The testing positivity $P$ can be expressed as: $$ P = \frac{Y T_Y}{T}=\frac{Y(1-B\Phi)}{T}=\frac{Y}{T}(1-\frac{(1-T)\Phi}{1-Y+Y\Phi})$$

[`HR_Test_positivity_figure.R`](HR_Test_positivity_figure.R): Generate figures for test positivity curves as a function of $\Phi=e^{-\phi}$ ($\phi$ is the hazard offset of infection), testing proportion and prevalence using Hazard Ratio idea.

-   [`./pix/HR_test_positivity_vs_phi-test_prop.png`](./pix/HR_test_positivity_vs_phi-test_prop.png): group by different test proportion

-   [`./pix/HR_test_positivity_vs_test_proportion-phi.png`](./pix/HR_test_positivity_vs_test_proportion-phi.png): group by different $\phi$ value

-   [`./pix/HR_heatmap.png`](./pix/HR_heatmap.png): Heat map used to show the boundary condition of $\Phi \geq 1-T/Y$.

### Alternative:

Treat $B$, $\Phi$ both as constant. Let $T$ as a function: $$T=T(B,\Phi,Y)=(1-Y)T_B+Y T_Y=(1-Y)(1-B)+Y(1-B\Phi)$$ Also, $P$ as a function: $$P=P(B,\Phi,Y)=\frac{Y T_Y}{T}=\frac{Y}{T}(1-B\Phi)=\frac{Y T_Y}{(1-Y)T_B+Y T_Y}=\frac{Y(1-B\Phi)}{(1-Y)(1-B)+Y(1-B\Phi)} $$

And try to use both time series of $T$ and $P$, trying to inference $(B,\Phi,Y(..))$ with MLE.

(??) Make $P$, $Y$ mutually independent in expression for fitting. $$T=(1-Y)(1-B)+Y(1-B\Phi) \Rightarrow Y=\frac{ 1-B - T }{B(1 - \Phi)}$$

$$P=\frac{Y}{T}(1-B\Phi)=\frac{(1-B-T)(1-B\Phi)}{TB(1-\Phi)} $$ Seems does not really matter in fitting.

## Odds Ratio

New updates from @dushoff on Nov 19th on slides

-   Hazard
    -   Pro(+): Usage is relatively easy to justify.
    -   Con(-): Probability curves characterized by hazard has steep drops near $0$ and "suddenly terminate" at $0$
-   Log Hazard
    -   Pro(+): To avoid con of hazard, log(Hazard) is used with a better look on neighborhood near $0$ and negative intervals.
    -   Con(-): It has unfavorable behavior in positive interval ($P \rightarrow 1$ very fast near $2$).
    -   Con(-): Harder to justify usage
-   Log Odds:
    -   Pro(+): Better behavior on both positive and negative interval
    -   Pro(+): Behave close to hazard in positive interval and close to log hazard in negative interval
    -   Pro(+): Still fairly easy to justify
    -   Pro(+): (??) More flexible/nonlinear behavior
    -   Con(-): Harder to compute

### Assumptions and notations

If we assume a constant odds ratio between test probabilities of positive and negative people, this should also parameterize a sensible set of curves.

We use similar notation with Hazard idea:

-   $T$ is the testing proportion of the whole population (also the average testing probability).

-   $P$ is the testing positivity.

-   $Y$ is the prevalence (also the average infection probability)

-   $B$ is the corresponding odds of uninfected/negative being tested.

    -   $T_B$ be the baseline probability of uninfected/negative individual being tested. $$B=\frac{T_B}{1-T_B} \Leftrightarrow T_B=\frac{B}{1+B}$$

    -   **Constraint** $B>0$

    -   (??) Maybe for later we use $b=log(B)$ as parameterization

    -   (??) Similar like hazard ratio, should $B$ be a constant as a result of behavior factor of testing? Or $B$ is a intermediate variable?

-   $\Phi$ is the (constant) odds ratio (strength of the association) between test probabilities of positive and negative people. $$\Phi=\frac{\frac{T_Y}{1-T_Y}}{\frac{T_B}{1-T_B}}=\frac{\frac{T_Y}{1-T_Y}}{B} \Leftrightarrow \Phi B= \frac{T_Y}{1-T_Y} \Leftrightarrow T_Y=\frac{B\Phi}{1+B\Phi}$$

    -   $T_Y$ is the probability of infected/positive individual being test.

    -   **Constraint** $\Phi>1$

    -   The probability $T_B$ is monotonically increase with corresponding odds.

    -   (??) $\phi=log(\Phi)$ later for better parameterization, $\phi>0$

### Target

Same with all other model: We would like to use time series of $T$ and $P$ as observed from data and trying to inference $Y$ and $\Phi$ with a statistical framework.

### Model

Similar with hazard model, we can express $T$ in: $$T = (1-Y) T_B + Y T_Y = (1-Y)\frac{B}{1+B}+Y\frac{B\Phi}{1+B\Phi}$$

We can use this expression to substitute $B$ as constant (or solve $B=B(T,Y,\Phi)$ as a variable). $$B=\frac{X \pm \sqrt{X^2+4T(1-T)\Phi}}{2(1-T)\Phi}$$ where $$X=(T-Y)\Phi+T+Y-1=+T(\Phi+1)-Y(\Phi-1)-1$$

-   Observation: no constraint on $X$, but behavior might change dramatically at $X=0$.
    -   Practically, we consider case that $T<Y, T+Y<1 \Rightarrow X<0$.
-   (TO DO??) Interpretation of $X$?

Only add(+sqrt) solution is feasible since $B>0$ as odds. Minus(-sqrt) solution is always negative.$$B=\frac{X + \sqrt{X^2+4T(1-T)\Phi}}{2(1-T)\Phi}$$ - (TO DO??) How $B$, thus $T_B$ change as a function of $T,Y,\Phi$?

Then the testing positivity $P$ can be expressed as: $$P=\frac{Y T_Y}{T}=\frac{Y}{T}\frac{B\Phi}{1+B\Phi}=\frac{Y}{T}\frac{X+\sqrt{X^2+4T(1-T)\Phi}}{X+\sqrt{X^2+4T(1-T)\Phi}+2(1-T)}$$ $$=\frac{Y}{T}(1+\frac{2(1-T)}{\sqrt{X^2+4T(1-T)\Phi}})^{-1}=\frac{Y}{T}(1-\frac{2(1-T)}{\sqrt{X^2+4T(1-T)\Phi}+2(1-T)})$$

[`or_mac.rmd`](or_mac.rmd): @dushoff description with preliminary calculation on maxima.

### Alternative:

Again, like hazard idea, treat $B$, $\Phi$ both as constant. Let $T$ as a function: $$T=T(B,\Phi,Y)=(1-Y)T_B+Y T_Y=(1-Y)(1-B)+Y(1-B\Phi)$$ Also, $P$ as a function: $$P=P(B,\Phi,Y)=\frac{Y T_Y}{T}=\frac{Y T_Y}{(1-Y)T_B+Y T_Y}=\frac{Y\frac{B\Phi}{1+B\Phi}}{(1-Y)\frac{B}{1+B}+Y\frac{B\Phi}{1+B\Phi}}=\frac{Y \Phi (1+B)}{1+(B+Y)\Phi-Y} $$

And try to fit both $T$ and $P$ as independent(???), trying to inference $(B,\Phi,Y(..))$ with MLE.

# Simulation and Fitting

## Simulation Phase:

### Exponential Growth Prevalence

Consider a population with size $N$. Starting with a exponential growth model for $Y$: $$Y(Y_0,r,t)=\frac{Y_0}{N} e^{r t} $$

-   $Y_0$ is the initial infection **number**
-   $r$ is the growth rate.

Assume we know the real value $B=\tilde{B}$, $\Phi=\tilde{\Phi}$, $\tilde{Y_0}$ and $\tilde{r}$, then we can construct time series of prevalence $\tilde{Y}$, expected testing positivity $\tilde{P}$ and expected testing proportion $\tilde{T}$.

We can then generate the data: time series for observed testing number $N T^{*}$ and observed number of positive tests $N T^{*} P^{*}$.

-   $T^{*}$ is the observed testing proportion

-   $P^{*}$ is the observed testing positivity

-   For observed number of test, we assume: $$ NT^{*} \sim \text{Pois}(N \tilde{T})$$

-   For observed number of positive test, we assume $$N T^{*} P^{*} \sim \text{Binom}(NT^{*},\tilde{P})$$

## Fitting Phase

Based on what we known:

-   Model equations

-   Data of $P^{*}$, $T^{*}$

-   Total population size

-   Poisson and Binomial distribution

Trying to find the best fit combination of ($\hat{B}$, $\hat{\Phi}$, $\hat{r}$ and $\hat{Y}_0$) using a MLE approach. $$ log(\mathcal{L}(\vec{T^{*}},\vec{P^{*}}|B,\Phi,r,Y_0))=\sum_{t=0} [log(\mathcal{L}_{Pois}(T^{*}(t)|NT(B,\Phi,Y(t;r,Y_0)))+log(\mathcal{L}_{binom}(P^{*}(t)|NT(B,\Phi,Y(t;r,Y_0)),P(B,\Phi,Y(t;r,Y_0)))]$$

## Related Positivity Bias Project with PHAC and PHO

[`simple.md`](simple.md) gives a brief introduction to current ideas of the positivity Bias project with David C. or Kevin B.
