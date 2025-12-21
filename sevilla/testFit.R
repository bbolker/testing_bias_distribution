library(shellpipes)
manageConflicts(add="dplyr")

library(purrr)
library(dplyr)
library(bbmle)

loadEnvironments()
obs <- rdsRead()

set.seed(seed)

mod <- mle2(testLike
	, start = list(beta=bstart, D=Dstart, I0=Istart, h=hstart)
	, data = list(obs=rdsRead()
		, steps=steps, deltat=deltat
	)
	, fixed = list(N=N)
)

summary(mod)
rdsSave(mod)
