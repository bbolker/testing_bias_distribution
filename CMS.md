
## Title
Improving infectious disease prevalence estimation and parameter inference using number of tests and positivity data

## Abstract 
Prevalence of infection is a critical variable for modeling infectious dynamics and public health decision-making. However, estimating true prevalence from surveillance data remains challenging. Here we discuss three novel attempts to model testing mechanism using the beta-distribution, hazard ratios and odds ratios respectively. These methods aim to link prevalence with the number of tests, test positivity, and test characteristics at each data point in a robust, flexible, and theoretically justified manner. We further present a data-fitting framework based on the odds ratio approach and demonstrate its performance using simulated datasets as a proof of concept.

## Introduction
Prevalence of the infection, i.e. proportion of infectious individuals in the population, is a critical state variable when modeling infectious dynamics.
For many model, to fit the parameters with data, we also to model relationship between prevalence as prediction and data as observation.  
From public health perspective, it is also the key information for understanding the status of the outbreak and for decision making.

However, it is hard to directly measure prevalence for large population or deduce it from the typical surveillance datasets.
In many surveillance datasets, we can only observe the portion of cases being tested and reported, but have little information about the general population.
One common approach assumes prevalence is equal to test positivity, which is the proportion of positive cases among all tested cases. 
Such assumptions implicitly indicates that the population being tested is a random sample of the general population.

But in real-world, the cohort of tested cases and untested population could be quite different, just due to the tests themselves is not random, which affect prevalence and positivity directly and significantly.

In general, infected people with symptoms is more likely to search medical care and be tested. 
And more specifically, the more severe the symptoms, the more likely someone is to be tested.
While for many infectious disease like influenza, sub-clinical cases might be dominant in the infected population, but might rarely being tested.

Testing will also depend on how they got infected, e.g. imported cases and close contacts of high profile infection like COVID, are more likely to be tested. 
On contrast, an asymptomatic, community-infected person could never be tested.

Therefore, we might want to consider more carefully why and how people are getting tested when estimating prevalence from surveillance data.
There are many factors might related to the testing and prevalence relation, like
- Delay among infection, being infectious and being tested
- Heterogeneity in testing behavior and test policies
- Impact of public awareness and medical
But here we focus on building basic, simple models for testing mechanic and priority on homogeneous population as a starting point.

The 3 models we discussed today, is focused on modelling test mechanics and priority with simple, uniform approaches characterized just by 1 to 2 parameters, that connect **prevalence** $Y$ to
- Number of tests done\Test proportion $T$
- Positive cases reported from the tests $P$
which are commonly collected into datasets by surveillance programs.


