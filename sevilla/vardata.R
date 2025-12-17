library(shellpipes)
library(purrr)
library(dplyr)

loadEnvironments()

set.seed(seed)

epi <- simulate(sirRates 
	, states = (list(t=0, S=N-I0, I=I0))
	, params = (list(beta=beta, D=D, N=N, deltat = deltat))
	, steps=100
)

summary(epi)

epi <- (epi
	|> mutate(t=t
		, V = 1 - (S+I)/N
		, a = w0 + wI*I*exp(-alpha*V)/N
		, rho = 1-exp(-a)
		, pos = rpois(nrow(epi), rho*I) 
		, neg = rpois(nrow(epi), rho*(N-I)/h)
		, pp = pos/(pos+neg)
	)
)

saveEnvironment()

## Fitting scripts only get obs
obs <- (epi |> select(t, pos, neg, pp))
summary(obs)
rdsSave(obs)

