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

## Calculate likelihood of parameters (beta, D, I0, h)
## for a data frame (obs) of pos/neg tests
## Don't use dpois for any of this, please!
testLike <- function(beta, D, I0, N, h, obs, steps, deltat, hmult=10){
	epi <- simulate(sirRates 
		, states = (list(t=0, S=N-I0, I=I0))
		, params = (list(beta=beta, D=D, N=N, deltat = deltat))
		, steps=steps
	)
	b <- calcBaseHazard(h, N, epi$I, obs$neg, obs$pos, hmult)
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
	## My probability of being a positive test is the product of 
	## prevalence and testing risk
	epi <- (epi
		|> mutate(
			, tprob = I*(1-exp(-h))/N
			, ll = dbinom(obs$pos, N, tprob, log=TRUE)
		)
	)
	## print(epi)
	return(-sum(epi$ll))
}

## Calculate the MLE floating baseline given relative hazard
## (i.e., return the base hazard that combines with provided rel to maximize the likelihood)
## rel is relative hazard
## N/I are population/incidence
## neg/pos are test results
## hmult provides a buffer for uniroot to bracket the solution
calcBaseHazard <- function(rel, N, I, neg, pos, hmult=10){
	opt <- function(x, rel, V, I, neg, pos){
		pN <- 1-exp(-x)
		pP <- 1-exp(-x-rel)
		return(
			(neg/pN - (V-neg)/(1-pN))*(1-pN)
			+ (pos/pP - (I-pos)/(1-pP))*(1-pP)
		)
	}
	hN <- -log(1-neg/V)/hmult
	hP <- -log(1-pos/I)*hmult
	u <- uniroot(opt, c(hN, hP-rel), rel=rel, V=N-I, P=I, neg=neg, pos=pos)
	return(u$root)
}

saveEnvironment()
