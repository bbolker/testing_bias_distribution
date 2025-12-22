## Just to test the math behind calcBase. Notation is outdated
library(bbmle)

binLike <- function(h, n, i){
	p <- 1-exp(-h)
	return(i*log(p) + (n-i)*log(1-p))
}

obsnll <- function(base, rel, N, P, neg, pos){
	return(
		- binLike(base, N, neg)
		- binLike(rel+base, P, pos)
	)
}

fitBase <- function(rel, N, P, neg, pos, start=1){
	fit <- mle2(obsnll
		, start = list(base=start)
		, fixed = list(rel=rel)
		, data = list(N=N, P=P, neg=neg, pos=pos)
	)
	b <- coef(fit)[["base"]]
	return(1 - exp(-b))
}

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

system.time(replicate(1000, fitBase(0, 30, 10, 3, 2)))
system.time(replicate(1000, calcBase(0, 30, 10, 3, 2)))
