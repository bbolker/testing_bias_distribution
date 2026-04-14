The following work is based on the hazard model with a time varying, free baseline hazard $b$ for uninfected population and a constant risk offset $\phi$ caused by infection.

For simplicity, we denote $B(t)=e^{-b(t)}$ and $H=e^{-\phi}$, s.t. negative testing probability is denote by $p_-=T_B=1-B(t)$ and positive testing probability is $p_+=T_Y=1-H\times B(t)$.

## Likelihood of $B$
See [spainMath.md](spainMath.md) for more details.

Let's say that group $k$ has size $n_k$ of whom $t_k$ are tested, and has hazard multiplier $\Phi$ (so that $\Phi_+ = H$ and $\Phi_- = 1$).

The probability of $t_k$ tests, given the parameters, is
$$ P_k = \binom{n_k}{t_k} (1-\phi_kB)^{t_k} (\phi_kB)^{n_k-t_k}$$
$$ \propto (1-\phi_kB)^{t_k} B^{n_k-t_k} \equiv L_k,$$
after dropping factors that don't depend on $B$.

The overall likelihood is:
$$\begin{aligned}
	L & \propto \prod_k{(1-\phi_kB)^{t_k} B^{n_k-t_k}}
	\\
	& =  \prod_k{(1-\phi_kB)^{t_k} B^{\sum(n_k-t_k)}}
	\\
	& =  \prod_k{(1-\phi_kB)^{t_k} B^{N-T}}
\end{aligned}$$
where $T$ is the total population of tested individuals.

Thus, the maximum-likelihood estimate of $B$ does not depend on the group sizes $n_k$, but only their sum the population size. 
This is because in this parameterization, an increase in $B$ makes each negative test more likely in the same proportion, thus has the same relative effect on their product no matter how they are distributed. 
Note that the $n_k$ **do** affect the original probability, via the $\prod \Phi_k^{n_k}$ that we scaled out.

We can then write log likelihood
$$ \ell = \sum_k t_k \log(1-\Phi_kB) + (N-T) \log(B)$$
and look for zeros of 
$$ \frac{d\ell}{dB}=\ell'_B = - \sum_k \phi_k \frac{t_k}{1-\phi_kB} + (N-T)/B$$

In terms of our original two groups, we could write this for example as:
$$ \ell'_B = - \frac{x}{1-B} - H\frac{y}{1-HB} + \frac{N-x-y}{B}=0,$$
where $x$; $y$ is the number of negative and positive tests observed.
$$\Leftrightarrow  N H B^2 - ((N - x)H + N - y)B +(N - x - y) =0$$
As a result, we have two potential $B$ solution that maximize the likelihood $$B=\frac{(N-x)H+N-y \pm \sqrt{((N-x)H+N-y)^2-4N H (N-x-y)} }{2NH}$$
As $(N-x)H+N-y=NH+N-xH-y > NH$ and $H\in(0,1)$, there is a good chance that the $B_{+}$ solution will be larger than $1$.


## Fitting Machine
### Fixed Baseline Test Cases
`OR_Sim_SIR.R`: mle2 fitting for basic SIR model with fixed baseline
`OR_Sim_SIR_MacPan.R`: macpan2 fitting for basic SIR model with fixed baseline
`OR_Sim_SIR_Season.R`: macpan2 fitting for seasonal SIR model with fixed baseline

All these machine is using the binomial version simulation and likelihood function:
#### For data
Prevalence based on infectious dynamics $Y = I/N$  
Expected test proportion $T = (1-Y)*T_B+Y*T_Y$
Expected test positivity $pos = Y \times T_Y/T$              
Observed test proportion $OT \sim \text{Binom}(N,T)$
Observed test positivity $OP \sim \text{Binom}(OT,pos)$

Only $OP(t), OT(t)$ is observed as data

#### For fitting
Prevalence based on infectious dynamics $\hat{Y} = I/N$  
Expected test proportion $\hat{T} = (1-\hat{Y})*\hat{T}_B+\hat{Y}*\hat{T}_Y$
Expected test positivity $\hat{pos} = \hat{Y} \times \hat{T}_Y/\hat{T}$
The fitting objective/likelihood function:
$$\ell = - \sum_{t} (\text{dbinom}(OT(t), N, \hat{T(t)})) - \sum_t(\text{dbinom}(OP, OT, \hat{pos}))$$

### Float Baseline Fitting
Current modularized macpan2 version
See `Makefile` Modularized float-baseline fitting 2026 Mar 08 section for more details of .R files and flow 

Flexible baseline setting:
$$b(t) = w_0+w_I \times I(t)/N \times e^{(-\alpha(1-S(t)/N))}$$
with $B= e^{-b}$, $\Phi=e^{-h}$, $T_B=1-B$ and $T_Y=1-B\Phi$
### Binomial Version
The simulation of data is basically the same idea with fixed baseline, but adding $b(t)$ to calculate $B$ and $T_B$.

For the simulator of fitting procedure, compare to fixed baseline version the following part is added/modified for floating baseline $b$:
- take observed $OP(t)$ and $OT(t)$ to calculate positive and negative case at time $t$ for calculating $B_{lik}$ and take all such result from data as fixed parameters.
	`pos ~ time_var(pos_values, pos_changepoints)`
	`neg ~ time_var(neg_values, neg_changepoints)`
Then calculate $B_{lik}$ based on previous likelihood derivation:
$B_{lik} = \frac{1}{2 N \Phi}\sqrt{((N-neg)\Phi+N-pos) - (((N-neg)\Phi+N-pos)^2-4N\Phi(N-pos-neg))}$
then treat $T_B = 1 - B_{lik}, T_Y = 1 - \Phi B_{lik}$
The rest is the same mechanism with fixed baseline approach.

The fitting result:
- `const.fixed.check.Rout.pdf` Constant baseline data with fixed baseline fitting, works well in general
- `const.flex.check.Rout.pdf` Constant baseline data with flexible baseline fitting, not working well. The negative log-likelihood is lower than fixed version, but the scale of prevalence is off with the "real" case in data simulation.
- `float.flex.check.Rout.pdf` Float baseline data with flexible baseline fitting, not working for the same off-scaling problem

### Poisson Version
#### JD's Comment: 
We are currently modeling the outbreak deterministically, meaning that we don't have an integer number of cases. Given that, maybe we should try using Poisson instead of binomial likelihoods. 




#### BB's Comment: 
My concern, which might be unfounded, is that we need both test positivity (+ tests/total tests) and testing fraction (tests/pop size) to be small in order for this to work.


RZ's conjecture would be the $B_{lik}$ is independent of prevalence $Y$ from the outbreak, given too much degree-of-freedom to the system. 

Even with constant $\Phi$, a prevalence $Y$ curve with different scale could still have the same or even better observed number of cases as the true $Y$ curve, since $B$ is very flexible and not related to scale of $Y$ at any single point. So the fitting result keeps the shape of $Y$ but mess up the scale of $Y$
