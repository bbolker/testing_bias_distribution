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

testLike <- function(beta, D, I0, N, h, obs, steps, deltat, hmult=10){
	epi <- simulate(sirRates 
		, states = (list(t=0, S=N-I0, I=I0))
		, params = (list(beta=beta, D=D, N=N, deltat = deltat))
		, steps=steps
	)
	b <- calcBase(h, N, P, neg, pos, hmult)
	## What????
	epi <- (epi
		|> mutate(
			, ll = dpois(obs$pos, posPred, log=TRUE)
		)
	)
	return(-sum(epi$ll))
}

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

## Calculate the MLE floating baseline given relative hazard
calcBase <- function(rel, N, P, neg, pos, hmult=10){
	opt <- function(x, rel, N, P, neg, pos){
		pN <- 1-exp(-x)
		pP <- 1-exp(-x-rel)
		return(
			(neg/pN - (N-neg)/(1-pN))*(1-pN)
			+ (pos/pP - (P-pos)/(1-pP))*(1-pP)
		)
	}
	hN <- -log(1-neg/N)/hmult
	hP <- -log(1-pos/P)*hmult
	u <- uniroot(opt, c(hN, hP-rel), rel=rel, N=N, P=P, neg=neg, pos=pos)
	invisible(1-exp(-u$root))
}

saveEnvironment()
