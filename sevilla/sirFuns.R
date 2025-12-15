library(shellpipes)

sirRates <- function(states, params){
	with(c(states, params), {
		inc <- beta*S*I/N
		rec <- I/D
		return(c(-inc, inc-rec))
	})
}

## keep states as a list, but accumulate it as a vector
simulate <- function(rateFun, states, params, steps) {
	v0 <- unlist(states, use.names = TRUE)
	p <- length(v0)
	out <- matrix(NA_real_, nrow = steps + 1, ncol = p)
	colnames(out) <- names(v0)
	out[1, ] <- v0

	for (t in 1:steps) {
		grad <- c(1, rateFun(states, params))
		v <- unlist(states, use.names = FALSE)
		vnew <- v + grad * params$deltat

		out[t + 1, ] <- vnew

		states <- as.list(vnew)
		names(states) <- colnames(out)
	}

	return(as.data.frame(out))
}

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
