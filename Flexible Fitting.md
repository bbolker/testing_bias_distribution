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
The fitting simulation is the same as 



### Float Baseline Fitting
Current modularized macpan2 version 
`parameter`
### Binomial Version
