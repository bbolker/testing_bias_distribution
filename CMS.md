
## Title
Improving Infectious Disease Prevalence Estimation and Parameter Inference Using Number of Tests and Positivity Data

## Abstract 
Prevalence of infection is a critical variable for modeling infectious dynamics and public health decision-making. However, estimating true prevalence from surveillance data remains challenging. Here we discuss three novel attempts to model testing mechanism using the beta-distribution, hazard rates and odds ratios respectively. These methods aim to link prevalence with the number of tests, test positivity, and some test characteristics at each data point in a robust, flexible, and theoretically justified manner. We further present a data-fitting framework based on the odds ratio approach and demonstrate its performance using simulated datasets as a proof of concept.

## Introduction
Prevalence of the infection $Y(t)$, i.e. proportion of infectious individuals in the population at time $t$, is a critical variable when modeling infectious dynamics, as it often determines the force of infection.

Also, when fitting models the parameters with data, we often need to consider relationship between prevalence as model prediction and data as observation.  
From public health perspective, it is also the key information for understanding the status of the outbreak and for decision making.

However, it is hard to directly measure prevalence for large population or deduce it from the typical surveillance datasets. 
A key reason for such difficulty is that, we usually are only able to observe cases being tested and reported, but have little information about the untested population.

One common approach assumes prevalence is equal to test positivity, which is the proportion of positive cases among all tested cases. 
Such assumptions implicitly indicates that the tests are distributed totally, uniformly random to the whole population.

But in real-world, most tests is not distributed at random, and how test is distributed should affect prevalence and positivity directly and significantly.

On the other hand, infected people with symptoms is more likely to search medical care and be tested. 
More specifically, the more severe the symptoms, the more likely someone is to be tested.
While for many infectious disease like influenza, sub-clinical cases might be dominant in the infected population, but might rarely being tested.

Testing will also depend on how they got infected, e.g. imported cases and close contacts of high profile infection like COVID, are more likely to be tested. 
On contrast, an asymptomatic, community-infected person could never be tested.

A common 

Therefore, we might want to consider more carefully why and how people are getting tested when estimating prevalence from surveillance data.
There are many factors might related to the testing and prevalence relation, like
- Delay among infection, being infectious and being tested
- Heterogeneity in testing behavior and test policies
- Impact of public awareness and medical
But here we focus on building basic, simple models for testing mechanic and priority on homogeneous population as a starting point.

The 3 models we discussed today, is focused on modelling test mechanics and priority with simple, uniform approaches characterized just by 1 to 2 parameters, that connect **prevalence** $Y$ to
- Proportion of tested individuals among whole population $T$ 
- The test positivity $P$
in a short time duration, which could be the reporting period of surveillance data.
If assume we have a large enough population $N$, then we can easily calculate such proportions from normalize the number of conducted tests and number of positive tests, which are commonly collected in surveillance data.

Also, we assume the total number of positive tests is small compare to the whole population, i.e. $0 < T\ll$ 1, thus the positive cases being tested has very limited impact to the infectious dynamics in the whole population.

### Beta-distribution Model
We start with a top-down approach, consider how tests are distributed in the whole population. This approach is intuitively more proper for the situation that test opportunities are limited.

There are two extreme scenarios that are easy to reason about:
- If testing is totally random, the proportion of positive tests should be the same as the prevalence. 
	- Increasing the number of tests when the prevalence is constant leads to a proportional increase in the number of cases. 
	- If the prevalence is increasing exponentially with rate $r_1$ and testing is increasing at rate $r_2$, the number of cases increases at rate $r_1+r_2$.
- If testing is perfectly focused on infected people, then tests are 100% positive until we run out of infected people (i.e. $T$ > $Y$). The proportion of positive tests is $\min(1,Y/T)$ .

Now let’s consider that in general we test people in (approximate) order of their probability of infection (ignoring tests skipped because of constraints). 
We can think about: 
- A distribution in the population of probability of infection
	- E.g. higher for known contacts of infected people, travelers from high-risk areas, etc. 
- A distribution in the probability of testing (which ideally lines up with the risk). 

We could think about having two distributions, but these can’t really be separated, so we tried to use a single beta-distribution to combine these two probability, with PDF: $$f_\beta(x;a,b)=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}x^{a-1}(1-x)^{b-1}=\frac{x^{a-1}(1-x)^{b-1}}{B(a,b)}$$
- where $B(a,b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}$ is the beta-function.
- $f_\beta$ gives the population density of probability $x$ for individuals being infected.
- we could parameterized the distribution by
	- Its mean $E[X]=\frac{a}{a+b}=Y$ equals to the prevalence.
	- Its dispersion $\phi=1-e^{-\frac{1}{a+b}}$ that characterize the "testing focus"
	- Such that $$a=\frac{-Y}{ln(1-\phi)}, b=\frac{1-Y}{ln(1-\phi)}$$
If we always test 'from the top down', i.e. calculate the mean value of the top fraction $T/N$ of the population distribution,
- As $\phi \rightarrow 1$ we end up with point masses $1-Y$ on 0 and $Y$ on 1, correspond to perfect-target testing (cases=$\min(T,Y)$)
- As $\phi \rightarrow 0$ we end up with a point mass on $Y$, correspond to purely random testing (cases=$\textrm{Binomial}(i,T)$). 

Example beta-curve $Y=0.25, \phi=0.01$![[beta_phi001.png]]
Example beta-curve $Y=0.25, \phi=0.95$![[beta_phi095.png]]
Then using the distribution, we could construct the relation between $Y$ as mean of the distribution and the expected testing positivity, as the mean density of the tested portion (grey area):$$\begin{align}
\bar{P}&=\frac{\left(\int_{Q_\beta(1-T)}^1 f_\beta(x;Y,\phi) x \, dx\right)}{1-F_\beta(Q_\beta(1-T))} 
\\
& = \frac{\left(\int_{Q_\beta(1-T)}^1 f_\beta(x;Y,\phi) x \, dx\right)}{T}
\\
& = Y(1+\frac{Q_\beta(1-T)^a(1-Q_\beta(1-T))^b}{T a\mathbb{B}(a,b)})
\end{align}$$
- $F_\beta$ is the corresponding CDF of $f_\beta(x;Y,\phi)$
- $Q_\beta(z)$ is the $z$-quantile of the distribution
- $\mathbb{B}(a,b)$ is the complete beta function

![[test_positivity_vs_test_proportion-phi.png]]

In `R`, the built in `qbeta` function has numerical inconsistency, such that as the actual value of $Q_\beta(1-T)\rightarrow 0$, the output overflows to $0$
- In the framework, this happens when $\phi \rightarrow 1$ or $Y \rightarrow 1$
- Area gets larger when $T$ is small (not idea as we consider limited tests)
- Bad for data fitting procedure
![[log-diff_est-log.png]]
An alternative numeric approach we have tried, is to use linear estimation as $Q_\beta \rightarrow 0$, this leads to $$\bar{P} \approx \frac{b+1-(1-Q_\beta(1-T))ab}{b+1-(1-Q_\beta(1-T))(1-a)b}\approx\frac{b+1-ab}{b+1-(1-a)b}$$
- Numerical estimation of $\bar{P}$ that reduce the impact of the overflow
- Still not ideal for fitting: numerical cusp when shifting expression

Pros and Cons
- ($+$) Testing process is characterized by one parameter $\phi$
- ($-$) Hard to interpret, measure and change $\phi$ (
	- Can only get from fitting
	- Not ideal for fitting a changing $\phi$ 
	- Hard to Modeling testing policy change
- ($-$) Numerical overflow 
- ($-$) Performance concern for fitting data to dynamics
	- Prevalence $Y$ changing will change the distribution
### Hazard Rate Model
### Odd Ratio Model