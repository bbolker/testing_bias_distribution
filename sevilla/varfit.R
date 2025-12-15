library(shellpipes)
manageConflicts(add="dplyr")
loadEnvironments()

library(purrr)
library(dplyr)
library(bbmle)
set.seed(seed)

loadEnvironments()
obs <- rdsRead()
summary(obs)

mod <- mle2(varLike
	, start = list(beta=0.3, D=6, I0=4, N=N, h=1)
	, data = list(pos = obs$pos, neg = obs$neg
		, steps=steps, deltat=deltat
	)
	, fixed = list(N=N)
)

summary(mod)
rdsSave(mod)
