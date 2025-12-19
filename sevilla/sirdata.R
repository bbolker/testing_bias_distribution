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

summary(epi)

obs <- (epi
	|> transmute(t=t
		, pos = rpois(nrow(epi), rho*I) 
		, neg = rpois(nrow(epi), sig*S)
	)
)

summary(obs)

rdsSave(obs)

