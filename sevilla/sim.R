library(shellpipes)
library(purrr)
library(dplyr)

loadEnvironments()

set.seed(seed)

epi <- simulate(sirRates 
	, states = (list(t=0, S=N-I0, I=I0))
	, params = (list(beta=beta, D=D, N=N, deltat = deltat))
	, steps=steps
)

## Concern goes up with prevalence, and down with cumulative illness
epi <- (epi 
	|> mutate(base=concernFun(S, I, N, w0, wI, alpha))
)

summary(epi)

obs <- testResults(epi, hazFun, binPick, pars=list(hr=h), N)

summary(obs)

rdsSave(obs)

saveEnvironment()
