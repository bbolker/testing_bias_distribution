library(shellpipes)

sirRates <- function(states, params){
	with(c(states, params), {
		inc <- beta*S*I/N
		rec <- I/D
		return(c(-inc, inc-rec))
	})
}

simulate <- function(rateFun, states, params, steps) {
	v0 <- unlist(states, use.names = TRUE)
	p <- length(v0)

	out <- matrix(NA_real_, nrow = steps + 1, ncol = p)
	colnames(out) <- names(v0)

	# store initial state
	out[1, ] <- v0

	# main loop
	for (t in 1:steps) {
		grad <- c(1, rateFun(states, params))
		v <- unlist(states, use.names = FALSE)
		vnew <- v + grad * params$deltat

		# write the row
		out[t + 1, ] <- vnew

		# push back into list form if you need it for rateFun next step
		# this keeps names aligned with colnames
		states <- as.list(vnew)
		names(states) <- colnames(out)
	}

	# convert to data.frame only once at the end
	df <- as.data.frame(out, check.names = FALSE)
	return(df)
}

saveEnvironment()
