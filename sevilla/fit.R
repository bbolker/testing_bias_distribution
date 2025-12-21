library(shellpipes)
manageConflicts(add="dplyr")

library(purrr)
library(dplyr)
library(bbmle)

loadEnvironments()
obs <- rdsRead()

set.seed(seed)

cn <- 0.2

bstart <- rlnorm(1, log(beta), cn)
Dstart <- rlnorm(1, log(D), cn)
Istart <- rlnorm(1, log(I0), cn)
hstart <- rlnorm(1, log(h), cn)

## Key fitting function should be testLike, remember to use some sort of cool math thing that you already found.
## Also introduce prevLike, where pos/neg assumed drawn from true prevalence?

## ???? I'm developing here on plane, but these functions I guess should be elsewhere

mod <- mle2(posLike
	, start = list(beta=bstart, D=Dstart, I0=Istart, h=hstart)
	, data = list(obs=rdsRead()
		, steps=steps, deltat=deltat
	)
	, fixed = list(N=N)
)

summary(mod)
rdsSave(mod)
