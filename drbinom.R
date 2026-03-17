
## This could be made more efficient with a purpose-built vectorization function; maybe if we are worried about performance at some point
drbinom0 <- function(x, size, prob, log = TRUE, addmult=10){
	## if(is.null(eps)) eps <- prob
	if (x<size) return(dbinom(x, size, prob, log))
	adjexp <- size+addmult*(x-size)
	if (log) return(adjexp*log(prob))
	return(prob^(adjexp))
}

drbinom <- Vectorize(drbinom0, c("x", "size", "prob"))

data.frame(
	fragile=dbinom(1:13, 8, 0.8, log=TRUE)
	, robust=drbinom(1:13, 8, 0.8)
)
