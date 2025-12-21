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

concernFun <- function(S, I, N, w0, wI, alpha){
	 x = I/N
	 return(w0 + wI*x*exp(-alpha*(1-x)))
}

## A risk fun for additive hazards
hazFun <- function(b, pars){
	with(pars, return(data.frame(
		S = 1-exp(-b)
		, I = 1-exp(-b-hr)
	)))
}

saveEnvironment()
