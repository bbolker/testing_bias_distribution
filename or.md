
The point of the beta truncation approach is simply to parameterize sensible curves that connect test positivity to prevalence. There are probably there are simpler ways to do this, one of which was suggested by Ben's thoughts about logistic models.

If we assume a constant odds ratio between test probabilities of positive and negative people, this should also parameterize a sensible set of curves. 

I am going to use P, T, and V for prevalence, proportion tested, and test positivity. I may also use p, τ, and v for the corresponding odds ratios. b and rb are the odds ratios of negative and positive people getting tested. Given V and r we can solve for the value of b that matches a desired value of T, and thus calculate P. We will then want to invert this relationship so that we could take an observed value of test positivity P, and use testing rate T and our parameter r to calculate prevalence V.

If we know b and r, we can calculate proportion tested as 

