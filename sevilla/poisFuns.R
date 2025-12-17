library(shellpipes)

## To do: have a shared component of the downstream calculations here
## and in <foo>data.R
## Or not -- would that only work for sir and not here?
varLike <- function(beta, D, I0, N, h, pos, neg, steps, deltat){
	epi <- simulate(sirRates 
		, states = (list(t=0, S=N-I0, I=I0))
		, params = (list(beta=beta, D=D, N=N, deltat = deltat))
		, steps=steps
	)
	epi <- (epi
		|> mutate(
			, poissonRatio = (pos+neg)/(h*I+(N-I))
			, posPred = h*I*poissonRatio
			, negPred = (N-I)*poissonRatio
			, llp = dpois(pos, posPred, log=TRUE)
			, lln = dpois(neg, negPred, log=TRUE)
		)
	)
	return(-sum(epi$llp)-sum(epi$lln))
}

sirLike <- function(beta, D, I0, N, rho, obs, steps, deltat){
	epi <- simulate(sirRates 
		, states = (list(t=0, S=N-I0, I=I0))
		, params = (list(beta=beta, D=D, N=N, deltat = deltat))
		, steps=steps
	)
	epi <- (epi
		|> mutate(
			, pred = rho*I
			, ll = dpois(obs, pred, log=TRUE)
		)
	)
	# print(epi)
	return(-sum(epi$ll))
}

saveEnvironment()
