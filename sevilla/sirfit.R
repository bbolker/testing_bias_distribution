library(shellpipes)
manageConflicts(add="dplyr")

library(purrr)
library(dplyr)
library(bbmle)

loadEnvironments()
obs <- rdsRead() |> pull(pos)

set.seed(seed)

mod <- mle2(sirLike
	, start = list(beta=0.3, D=6, I0=4, N=N, rho=0.1)
	, data = list(obs= rdsRead() |> pull(pos)
		, steps=steps, deltat=deltat
	)
	, fixed = list(N=N)
)

summary(mod)
rdsSave(mod)
