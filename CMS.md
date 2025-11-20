
## Title
Improving infectious disease prevalence estimation and parameter inference using number of tests and positivity data

## Abstract 
Prevalence of infection is a critical variable for modeling infectious dynamics and public health decision-making. However, estimating true prevalence from surveillance data remains challenging. Here we discuss three novel attempts to model testing mechanism using the beta-distribution, hazard ratios and odds ratios respectively. These methods aim to link prevalence with the number of tests, test positivity, and test characteristics at each data point in a robust, flexible, and theoretically justified manner. We further present a data-fitting framework based on the odds ratio approach and demonstrate its performance using simulated datasets as a proof of concept.

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
- 

### Beta-distribution Model
We start with a top-down approach, consider how tests are distributed in the whole population.

There are two extreme scenarios that are easy to reason about:
- If testing is totally random, the proportion of positive tests should be the same as the prevalence (or incidence - I’m not being precise about whether we are measuring a flow or a stock, but I don’t think it makes a difference), regardless of either. Increasing the number of tests when the prevalence is constant leads to a proportional increase in the number of cases. The proportion of positive tests is i and the number of cases (i.e., positive tests) is iT. (I’m using i, t for proportion of the population infected or tested and I, T as the number infected or tested: I=iN, T=tN.) If the prevalence is increasing exponentially with rate r1 and testing is increasing at rate r2, the number of cases increases at rate r1+r2 (and doubling time log(2)/(r1+r2)).
- if testing is perfectly focused on infected people, then tests are 100% positive until we run out of infected people (i.e. proportion of pop tested > prevalence). The proportion of positive tests is min(1,i/t) and the number of cases is min(T,I).

Now let’s consider that in general we test people in (approximate) order of their probability of infection (ignoring tests skipped because of constraints). We can think about a distribution in the population of probability of infection (higher for known contacts of infected people, travelers from high-risk areas, etc.), and a distribution in the probability of testing (which ideally lines up with the risk). We could think about having two distributions, but these can’t really be separated (FIXME: explain better!).
