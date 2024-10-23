## brute-force beta quantile (inverse CDF) function
## params as in qbeta()
Qbeta0 <- function(p,shape1,shape2,lower.tail=FALSE) {
  fn <- function(x) {pbeta(x,shape1,shape2,lower.tail=lower.tail)-p}
  uniroot(fn,interval=c(0,1))$root
}
## vectorized version 
Qbeta <- Vectorize(Qbeta0,c("p","shape1","shape2"))

## qbeta() works sometimes, Qbeta() works sometimes ...
## try qbeta(), switch to Qbeta() if it throws a warning
## vectorized manually because we need to detect warnings at
## the individual qbeta() evaluation level
Qbeta2 <- function(p,shape1,shape2,lower.tail=FALSE) {
  n <- max(length(p),length(shape1),length(shape2))
  p <- rep(p,length.out=n)
  shape1 <- rep(shape1,length.out=n)
  shape2 <- rep(shape2,length.out=n)
  res <- rep(NA,n)
  for (i in seq(n)) {
    res[i] <- tryCatch(qbeta(p[i],shape1[i],shape2[i],lower.tail=lower.tail),
                       warning=function(w) {
                         Qbeta0(p[i],shape1[i],shape2[i],lower.tail=lower.tail)
                       })
  }
  res
}

##' @param i prevalence
##' @param t testing proportion
##' @param phi dispersion
##' @param numint use numerical integration?
##' @param qfun beta quantile function to use
##' @param plot.it illustrate?
##' @param phiscale ?? # this means if we would like to re-scale gamma in (0,+Inf) to phi in (0,1)
## FIXME: remove vestigial (qfun, numint) options
## FIXME: test phi=0; add special case if necessary (and maybe phi=1?)
prop_pos_test0 <- function(i,t,phi,
                           numint=FALSE,qfun=Qbeta2,
                           plot.it=FALSE,
                           phiscale="constrained",
                           debug=FALSE
                           ) {
    if (debug) cat("p1: ",i,t,phi,"\n")  
    ## cat(phi,"\n")
    if (phiscale=="constrained") {
        ## transform from (0,inf) to (0,1)
        ## y = inverse(1-exp(-x)) -> -log(1-y)
       phi_0 <- phi
       phi <- -log(1-phi)
       if (is.nan(phi)) {
        cat(phi_0)
        }
    }
  
    i<- min(i,1) # treat i>1 as i=1 (exp growth of I makes i>1? or is it just the optimizer doing its thing?)
    
    a <- i/phi; b <- (1-i)/phi
    if (plot.it) curve(dbeta(x,a,b),from=0,to=1)
    ## need special logic for extreme phi values?
    if (debug) cat("p2: ",t,a,(1-i)/phi,"\n")
    lwr <- qfun(t,a,b,lower.tail=FALSE)
    if (debug) cat("lwr: ",lwr,"\n")
    if (plot.it) abline(v=lwr,col=2,lty=2)
    ## replace with pbeta(...,lower.tail=FALSE) with appropriate multiplier?
    val <- if (numint) {
               fn <- function(x) x*dbeta(x,a,b)
               integrate(fn, lower=lwr, upper=1)$value
           } else {
               pbeta(lwr,a+1,b,lower.tail=FALSE)*(a/(a+b))
           }
    if (debug) cat("val: ",lwr,"\n")
    return(val/pbeta(lwr,a,b,lower.tail=FALSE))
}
prop_pos_test <- Vectorize(prop_pos_test0,c("i","t","phi"))


## 'pure' version of prop_pos_test0, for RTMB use
prop_pos_test1 <- function(i,t,phi) {
    ## transform phi??
    a <- i/phi; b <- (1-i)/phi
    lwr <- qbeta(1-t,a,b)
    ## note RTMB doesn't have lower.tail = FALSE option
    val <- (1.0 - pbeta(lwr,a+1,b))*(a/(a+b))
    return(val/t)
}

prop_pos_test2 <- function(i,t,phi,
                           method="log",
                           qfun=Qbeta2,
                           plot.it=FALSE,
                           phiscale="constrained",
                           debug=FALSE
) {
  if (debug) cat("p1: ",i,t,phi,"\n")  
  ## cat(phi,"\n")
  if (phiscale=="constrained") {
    ## transform from (0,inf) to (0,1)
    ## y = inverse(1-exp(-x)) -> -log(1-y)
    phi_0 <- phi
    phi <- -log(1-phi)
    if (is.nan(phi)) {
      cat(phi_0,"\n")
    }
  }
  
  i<- min(i,1) # treat i>1 as i=1 (exp growth of I makes i>1? or is it just the optimizer doing its thing?)
  
  a <- i/phi; b <- (1-i)/phi
  if (plot.it) curve(dbeta(x,a,b),from=0,to=1)
  ## need special logic for extreme phi values?
  if (debug) cat("p2: ",t,a,(1-i)/phi,"\n")
  lwr <- qfun(t,a,b,lower.tail=FALSE)
  if (debug) cat("lwr: ",lwr,"\n")
  if (plot.it) abline(v=lwr,col=2,lty=2)
  ## replace with pbeta(...,lower.tail=FALSE) with appropriate multiplier?
  ### Options for "method"
  ### "int" (brute-force integration)
  ### "cdf" (pbeta version)
  ### "simp" (simplified explicit version)
  ### "log" (test: log trick for exponential in "simp")
  val <- switch(method,
                int = {
                    fn <- function(x) x*dbeta(x,a,b)
                    integrate(fn, lower=lwr, upper=1)$value
                },
                cdf = pbeta(lwr,a+1,b,lower.tail=FALSE)*(a/(a+b)),
                simp = a/(a+b)*(t+(lwr^a*(1-lwr)^b)/(beta(a,b)*a)),
                log = a/(a+b)*(t+(exp(a*log(lwr)+b*log(1-lwr)-log(beta(a,b)*a)))),
                est = t*(b+1-(1-lwr)*a*b)/(b+1-(1-lwr)*(a-1)*b)
                )
  if (debug) cat("val: ",val,"\n")
  return(val/t)  ## CHECK denominator!
}
prop_pos_test_new <- Vectorize(prop_pos_test2,c("i","t","phi"))

## prop_pos_test_new(i=0.8, t= 0.0012, phi = c(0.01, 0.99))
