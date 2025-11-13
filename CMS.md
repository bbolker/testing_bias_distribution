## Title
Improving infectious disease prevalence estimation and parameter inference using number of tests and positivity data

## Abstract 
Prevalence of infection is a critical variable for modeling infectious dynamics and public health decision-making. However, estimating true prevalence from surveillance data remains challenging. Commonly used assumption—that prevalence is a function with just positivity as variable—are questionable for some diseases, because tested individuals might not represent the general population and the testing mechanism could also change with time. Here we discuss three novel attempts to model testing mechanism using beta-distribution, hazard ratio and odds ratio respectively. These methods aim to link prevalence with the number of tests, test positivity, and testing characteristics at each data point in a more robust, flexible, and theoretically justified manner. We further present a data-fitting framework based on the odds ratio approach, demonstrate its performance using simulated datasets as a proof of concept.

## Introduction
Prevalence of the infection, i.e. proportion of infectious individuals in the population, is a critical state variable when modeling infectious dynamics.
For many model, to fit the parameters with data, we also to model mechanism between prevalence as prediction and data as observation.  
From public health perspective, it is also the key information for understanding the status of the outbreak and for decision making.

However, it is hard to directly measure prevalence for large population or deduce it from the typical surveillance data.

In many surveillance data, we can only observe the portion of cases being tested and reported, but have little information about the general population.
Common approaches for estimating prevalence is just assume it is equal to, proportional to test positivity or it is a function depends on positivity, which is the proportion of positive, infectious cases among all tested cases.
However, the difference between tested and general population could be quite different, especially for those diseases that infected individuals might only have mild or no symptoms while still being infectious.

Therefore, we might want to consider more factors than just positivity when modeling the testing process and estimating prevalence.


There are other perspective related to the testing and prevalence problems, like
- Delay among infection, being infectious and being tested
- Heterogeneity in testing behavior and policy
- Other data source
But here we focus on construct a 



A commonly used approach is treated