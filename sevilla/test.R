
library(shellpipes)
library(purrr)

loadEnvironments()

test <- simulate(sirRates 
	, states = (list(t=0, S=0.99, I=0.1))
	, params = (list(beta=0.2, D=10, N=1, deltat = 1))
	, steps=100
)

summary(test)

rdsSave(test)

