library(shellpipes)

## epi is a simulation frame with a baseline concern measure
## riskFun converts baseline and a difference parameter to probability
## pickFun selects individuals from a population, with mean of demand
## pars can carry parameters for either function
## N is the population size S and I are compared to
testResults <- function(epi, riskFun, pickFun, pars, N){
	vprob <- with(epi, riskFun(base, pars))
	negDemand <- vprob$S*epi$S
	posDemand <- vprob$I*epi$I
	testDemand <- posDemand + negDemand
	posProb <- posDemand/testDemand
	tests <- pickFun(N, testDemand/N)
	return(tibble(t = epi$t
		, pos=pickFun(tests, posProb)
		, neg=tests-pos
	))
}

## A binomial pick function; does not need any parameters, but they are passed for some sort of possible future generality.
## Seems silly?
binPick <- function(N, p, pars){
	return(rbinom(length(p), N, p))
}

## HERE: figure out how to calculate likelihood
## I guess it's all right except we need to pass the right rho?
## somethin like invHaz(b+h)
posLike <- function(beta, D, I0, N, h, obs, steps, deltat){
	epi <- simulate(sirRates 
		, states = (list(t=0, S=N-I0, I=I0))
		, params = (list(beta=beta, D=D, N=N, deltat = deltat))
		, steps=steps
	)
	rho <- 1-exp(-h)
	epi <- (epi
		|> mutate(
			, pred = rho*I
			, ll = dpois(obs$pos, pred, log=TRUE)
		)
	)
	# print(epi)
	return(-sum(epi$ll))
}


saveEnvironment()
