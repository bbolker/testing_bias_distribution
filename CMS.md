
## Title
Improving infectious disease prevalence estimation and parameter inference using number of tests and positivity data

## Abstract 
Prevalence of infection is a critical variable for modeling infectious dynamics and public health decision-making. However, estimating true prevalence from surveillance data remains challenging. Here we discuss three novel attempts to model testing mechanism using the beta-distribution, hazard ratios and odds ratios respectively. These methods aim to link prevalence with the number of tests, test positivity, and test characteristics at each data point in a robust, flexible, and theoretically justified manner. We further present a data-fitting framework based on the odds ratio approach and demonstrate its performance using simulated datasets as a proof of concept.

## Introduction
Prevalence of the infection, i.e. proportion of infectious individuals in the population, is a critical state variable when modeling infectious dynamics.
For many model, to fit the parameters with data, we also to model relationship between prevalence as prediction and data as observation.  
From public health perspective, it is also the key information for understanding the status of the outbreak and for decision making.

However, it is hard to directly measure prevalence for large population or deduce it from the typical surveillance data.

In many surveillance data, we can only observe the portion of cases being tested and reported, but have little information about the general population.
Common approaches assume prevalence is equal to test positivity, which is the proportion of positive among all tested cases, or it is a linear function of positivity.
Such assumptions implicitly indicates that the population being tested is a random sample of the general population.
However, the difference between tested and general population could be quite different, especially for diseases that most infected individuals might only have mild or no clear symptoms while still being infectious; or the symptom is similar to, like influenza.
For those diseases, people with more severe or clear symptoms are often more likely to be tested, 



Therefore, we might want to consider more carefully why and how people are getting tested when estimating prevalence from surveillance data.
There are many factors might related to the testing and prevalence relation, like
- Delay among infection, being infectious and being tested
- Heterogeneity in testing behavior and policy
But here we start with building basic models for testing mechanic with limited complexity, 
The 3 models we discussed today, 

