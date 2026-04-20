The following work is based on the hazard model with a time varying, free baseline hazard $b$ for uninfected population and a constant risk offset $\phi$ caused by infection.

For simplicity, we denote $B(t)=e^{-b(t)}$ and $H=e^{-\phi}$, s.t. negative testing probability is denote by $p_-=T_B=1-B(t)$ and positive testing probability is $p_+=T_Y=1-H\times B(t)$.

# Likelihood of $B$
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


# Fitting Machine
## 1. Binomial Version
## 1.1 Fixed Baseline Test Cases
`OR_Sim_SIR.R`: mle2 fitting for basic SIR model with fixed baseline
`OR_Sim_SIR_MacPan.R`: macpan2 fitting for basic SIR model with fixed baseline
`OR_Sim_SIR_Season.R`: macpan2 fitting for seasonal SIR model with fixed baseline
`const.fixed.check.Rout`: modularized version of Macpan2 fitting. 

All these machine is using the binomial version simulation and likelihood function:
#### 1.1.1 For data
```R
  pY ~ I/N                          ## Prevalence based on SIR
, pSus ~ S/N                        ## Susceptible proportion  
, b ~ w0+wI*pY*exp(-alpha*(1-pSus)) ## Concern floating baseline hazard
, T_B ~ 1-exp(-b)
, T_Y ~ 1-exp(-b-h)
, T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
, T_pos ~ pY*T_Y/T_prop             ## Expected test positivity
, OT ~ rbinom(N,T_prop)             ## Observed test proportion
, OP ~ rbinom(OT,T_pos)             ## Observed test positivity
```
Only $OP(t), OT(t)$ is observed as data

#### 1.1.2 For fitting
```R
  pY ~ I/N                          ## Prevalence based on SIR
, pSus ~ S/N                        ## Susceptible proportion  
, b ~ w0+wI*pY*exp(-alpha*(1-pSus)) ## Concern floating baseline hazard
, T_B ~ 1-exp(-b)
, T_Y ~ 1-exp(-b-h)
, T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
, T_pos ~ pY*T_Y/T_prop             ## Expected test positivity
```
The fitting objective/likelihood function:
$$\ell = - \sum_{t} (\text{dbinom}(OT(t), N, \hat{T}(t))) - \sum_t(\text{dbinom}(OP, OT, \hat{pos}(t)))$$
- $\hat{T}$ and $\hat{pos}$ is the expected value from fitting.
## 1.2 Float Baseline Fitting
For more details about current modularized macpan2 version, see `Makefile`: 
- `modularized float-baseline fitting 2026 Mar 08`
- **The Binomial version is currently commented out !!**

Flexible baseline setting for test cases:
$$b(t) = w_0+w_I \times I(t)/N \times e^{(-\alpha(1-S(t)/N))}$$
with $B= e^{-b}$, $\Phi=e^{-h}$, $T_B=1-B$ and $T_Y=1-B\Phi$

#### 1.2.1 For data
The simulation of data is basically the same idea with fixed baseline, but adding $b(t)$ to calculate $B$ and $T_B$.

#### 1.2.2 For fitting
For the simulator of fitting procedure, compare to fixed baseline version the following part is added/modified for floating baseline $b$:
- take observed $OP(t)$ and $OT(t)$ to calculate positive and negative case at time $t$ for calculating $B_{lik}$ and take all such result from data as fixed parameters.
```R
pos ~ time_var(pos_values, pos_changepoints)
neg ~ time_var(neg_values, neg_changepoints)
```
Then calculate $B_{lik}$ based on previous likelihood derivation:
$B_{lik} = \frac{1}{2 N \Phi}\sqrt{((N-neg)\Phi+N-pos) - (((N-neg)\Phi+N-pos)^2-4N\Phi(N-pos-neg))}$
then treat $T_B = 1 - B_{lik}, T_Y = 1 - \Phi B_{lik}$
The rest is the same mechanism with fixed baseline approach, include the likelihood/objective function for MLE

### 1.3 The fitting result:
- `const.fixed.check.Rout.pdf` Constant baseline data with fixed baseline fitting, works well in general
- `const.flex.check.Rout.pdf` Constant baseline data with flexible baseline fitting, not working well. The negative log-likelihood is lower than fixed version, but the scale of prevalence is off with the "real" case in data simulation.
- `float.flex.check.Rout.pdf` Float baseline data with flexible baseline fitting, not working for the same off-scaling problem

#### 1.3.1 RZ's conjecture for the problem
The $B_{lik}$ expression is independent of prevalence $Y$ from the outbreak, given too much degree-of-freedom to the system. 

Even with constant $\Phi$, a prevalence $Y$ curve with different scale could still have the same or even better-fit observed number of cases as the true $Y$ curve, since $B$ is very flexible and not related to scale of $Y$ at any single point. So the fitting result keeps the shape of $Y$ but mess up the scale of $Y$.

# 2 Poisson Version
#### JD's Comment: 
We are currently modeling the outbreak deterministically, meaning that we don't have an integer number of cases. Given that, maybe we should try using Poisson instead of binomial likelihoods. 
The current hybrid-binomial approach is not perfect because it does not handle depletion properly, nor does it correlate properly between how many are picked and who is picked. 
Of course, the Poisson has (at least) the first of these problems.  
  
Our concrete example is if we have very few infected individuals with  
very high testing probability. The hybrid approach can greatly over-estimate variation in how many get tested. It will also miss an emergent negative correlation between how many are tested and positivity (since the most likely unusual events will be driven by high or low numbers of susceptible individuals testing).

My other suggestion is that we do a Poisson version for now (using Poisson for the “target truth” as well as the fitting) as an easier approach with fewer things to worry about. If that works well in some regime, we can try to translate to the #2 approach in Ben's Comment (noting for the record that the problematic non-integer is the binomial denominator, not the prediction).

#### Ben's Comment: 
My concern, which might be unfounded, is that we need both test positivity (+ tests/total tests) and testing fraction (tests/pop size) to be small in order for this to work.

The hybrid-binomial approach:
1) (obs tests/pop size), (pos tests/obs tests)
is incorrect as pointed out by JD because of non-independence, but are true binomials based on observed, discrete counts.
2) (obs tests/pop size), (pos tests/*predicted* tests):  
	- problematic when observed pos tests > predicted tests
	- predicted tests are non-integer 
BB don't think this is actually a big problem. Replacing the factorials in the normalization constant with gammas is no big deal, although can't quickly find a theoretical justification of what this actually means ...

I wonder if some of these problems would be alleviated (especially for #2, which seems more correct) if we were fitting simulations in a more realistic regime (i.e. we are probably, ultimately, interested in situations with low overall incidence/prevalence and testing fraction... ??)  
  
`2)`is the goal, including the gammas. 

### 2.1 Modularized Fitting
For time limitation, we only test this machine in one case within the reasonable regime to verify if alternating hybrid-binomial to Poisson will fix the fitting problem.
- Parameter: `float.params.R` $N=100,000$, with a time-varying $b(t)$ with low magnitude and small amplitude  
**TODO: Verifying the machine in "unrealistic" regime as Ben pointed out**

For more details about current modularized macpan2 version, see `Makefile`: 
- `modularized float-baseline fitting 2026 Mar 08`
- Poisson version files will have `.pois` in name
#### 2.1.1 For data
Now the observation is positive case OPos and negative case ONeg instead of OP and OT
```R
, ONeg ~ rpois(T_B*(N-I))
, OPos ~ rpois(T_Y*I)
```

#### 2.1.2 For fitting
Still, take OPos and ONeg from data as parameters:
```R
, pos ~ time_var(pos_values, pos_changepoints)
, neg ~ time_var(neg_values, neg_changepoints)
```

For a flexible fitting, use $B_{lik}$ as estimation of $B$
```R
, B_lik ~ 1/(2*N*Phi)*(((N-neg)*Phi+N-pos) - (((N-neg)*Phi+N-pos)^2-4*N*Phi*(N-pos-neg))^(1/2))
#, B_lik ~ exp(-w0)     ## for fixed fitting
, T_B ~ 1 - B_lik
, T_Y ~ 1 - Phi*B_lik
, ONeg ~ T_B*(N-I)
, OPos ~ T_Y*I
```
Here `OPos, ONeg` is for fitting simulation, automatically renamed as the `Sim_OPos` and `Sim_ONeg` following Macpan2/RTMB naming scheme in the likelihood: 
```R
calibrator$simulator$replace$obj_fn(~ - sum(dpois(obs_OPos, sim_OPos)) 
                                    - sum(dpois(obs_ONeg, sim_ONeg))
```
Mathematically, the likelihood now is 
$$\ell_\text{pois}=-\sum_{t}dpois(OPos, \hat{OPos})-\sum_{t}dpois(ONeg, \hat{ONeg})$$
### 2.2 Fitting result
- `const.fixed.pois.check.Rout.pdf` Constant baseline data with fixed baseline fitting, works well in general
- `const.flex.pois.check.Rout.pdf` Constant baseline data with flexible baseline fitting, not working well. The negative log-likelihood is lower than fixed version, but the scale of prevalence is off with the "real" case in data simulation.
- `float.flex.pois.check.Rout.pdf` Float baseline data with flexible baseline fitting, not working for the same off-scaling problem

Unfortunately, the result is not changing qualitatively.

### Next Steps?
Try the stepwise, spline and smooth functional form flexible baseline.