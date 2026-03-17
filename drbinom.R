
drbinom0 <- function(x, size, prob, log = TRUE){
	if (x<size) return(dbinom(x, size, prob, log))
	return(dbinom(size, size, prob^(x/size), log))
}

drbinom <- Vectorize(drbinom0, c("x", "size", "prob"))

data.frame(
	fragile=dbinom(1:13, 10, 0.8, log=TRUE)
	, robust=drbinom(1:13, 10, 0.8)
)
