
## This could be made more efficient with a purpose-built vectorization function; maybe if we are worried about performance at some point
drbinom0 <- function(x, size, prob, log = TRUE, eps=1e-3){
	if(is.null(eps)) eps <- prob
	if (x<=size) return(dbinom(x, size, prob, log))
	if (log) return(size*log(prob)+(x-size)*log(eps))
	return(prob^size*eps^(x-size))
}

drbinom <- Vectorize(drbinom0, c("x", "size", "prob"))

data.frame(
	fragile=dbinom(1:13, 10, 0.8, log=TRUE)
	, robust=drbinom(1:13, 10, 0.8)
)

data.frame(
	fragile=dbinom(1:13, 10, 0.8)
	, robust=drbinom(1:13, 10, 0.8, log=FALSE)
)
