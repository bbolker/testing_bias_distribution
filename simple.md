---
title: "simple/null modeling of multiplex testing"
---

We should explore a very simple world where we specify prevalence of two pathogens that share a duplex test and show patterns of how changes in prevalence in pathogen A affect positivity totals and proportions in pathogen B.

Suppose we have prevalences/incidences for
two diseases, 
$P_A$ and $P_B$.

The baseline hazard for testing is $\beta_0$ (i.e, prob (testing/
uninfected) is $1 - \exp (-\beta_0)$). Being infected with disease $A$ or \beta
increases hazard by $\beta_A$ or $\beta_B$.
So the fraction of the population testing positive for $A$ is

$P_A (1-\exp (-(P_0 + P_A)))$ and the test positivity rate for $A$ is
$$
\frac{P_A (1 - \exp(-(\beta_0 + \beta_A)))}{
\left[(1-\exp (-\beta_0)) (1-PA-P_B) + (1 -\exp (-\beta_0 + \beta_A))
p_A + (1-\exp(-\beta_0 + \beta_B)) P_B \right]
}
$$

Next steps:

- incorporate test sensitivity and specificity in the expressions above
- get some rough data from David C. or Kevin B. if possible (i.e. basic time series of test positivity and number of positive tests by day or week for influenza, RSV, COVID)? (Realizing that we don't know much about the mixture of single/multiplex tests ...)
- back-of-the-envelope calculations of the expected effects of disease $B$ on positivity of $A$ (for, say, influenza/COVID) in this null-model case
- make $\beta_0$ a function of $P_A$ and $P_B$? (Heightened awareness, contact tracing, etc.?) 
- What about age dependence?
- What would we need to get crude estimates of the $\beta$ values??
