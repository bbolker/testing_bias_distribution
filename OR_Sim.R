
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)

# Set Seeds
set.seed(13519)

## Initial true values:
T_B <- 0.04               ## uninf testing prob
T_Y <- 0.5                ## inf testing prob
B <- T_B/(1-T_B)          ## baseline odds of testing
Phi <- (T_Y/(1-T_Y))/B    ## inf vs uninf testing odds ratio

print(B)
print(Phi)

Y_0 <- 1e-4      ## initial prevalence
N <- 1e6         ## pop size
NY_0 <- N*Y_0    ## initial number infected

r <- log(2)/3    ## growth rate (doubling time = 3)
tmax <- 59       ## max simulation time
t <- c(0:tmax)
pts <- length(t) ## number of time points

true_param <- c("log_B"=log(B),"log_Phi"=log(Phi),"logY_0"=log(Y_0),"r"=r)


## Simulate the data
dat <- tibble(t=t
	## , pY = pmin(Y_0*exp(r*t), 1) ## Exponential growth
	, pY = 1/(1+(1/Y_0-1)*exp(-r*t))
	## , NY = rbinom(pts, N, pY)
	, NY = round(N*pY)
	, posTests = rbinom(pts, NY, T_Y)
	, negTests = rbinom(pts, N-NY, T_B)
)
dat

matplot(dat$t, dat[,-1], type = "l", log = "y")
legend("center", col = 1:4, lty = 1:4,
       legend = names(dat)[-1])

long_dat <- (dat
	|> select(-pY)
	|> pivot_longer(-t)
)

print(ggplot(long_dat)
 	+ aes(t, value, color=name)
 	+ geom_line()
 	+ scale_y_log10()
)

### function to calculate negative log-likelihood:
LL <- function(log_B, log_Phi, logY_0, r, dat, N, tmax, debug = TRUE,
               debug_plot = FALSE, plot_sleep = 1) {
    Y_0 <- exp(logY_0)
    B <- exp(log_B)
    Phi <- exp(log_Phi)
    T_B <- B/(1+B)
    T_Y <- B*Phi/(1+B*Phi)
    t <- c(0:tmax)
    pts <- length(t)
    ## simulated time series
    sim <- tibble(t=t
                ##, pY = pmin(Y_0*exp(r*t), 1)
                , pY = 1/(1+(1/Y_0-1)*exp(-r*t))
                ### round here???
                , NY = round(N*pY)
                # , NY = rbinom(pts, N, pY)
                )
    if(max(sim$pY) == 1 || any(sim$NY<dat$posTests) || any((N-sim$NY)<dat$negTests) || any(N<sim$NY)) return(NA)

  if (any(sim$NY<dat$posTests)) {
      cat("Underestimated infected population, pos tests > infected population", "\n")
  }
  if (any((N-sim$NY)<dat$negTests)) {
      cat("Overestimated infected population, neg tests > uninfected population", "\n")
  }
  postest_nll <- -sum(dbinom(dat$posTests, sim$NY, T_Y, log = TRUE))
  negtest_nll <- -sum(dbinom(dat$negTests, N-sim$NY, T_B, log = TRUE))
  out <- postest_nll + negtest_nll
  if (debug) {
      cat(B, Phi, logY_0, r, postest_nll, negtest_nll,
          out, "\n")
  }
  if (debug_plot) {
      par(mfrow= c(1,2), las = 1)
      ylim <- range(c(dat$posTests, dat$negTests,
                             sim$NY*T_Y, (N-sim$NY)*T_B))
      matplot(dat$t, dat[c("posTests", "negTests")], type = "p",
              pch = 1:2, log = "y",
              ylim = ylim)
      matlines(dat$t, cbind(sim$NY*T_Y, (N-sim$NY)*T_B))
      LLhist <<- c(LLhist, out)
      plot(LLhist - min(LLhist) + 1e-3, type = "b", log = "y")
      Sys.sleep(plot_sleep)
  }
  return(out)
}

real_ML <- LL(log(B),log(Phi),log(Y_0),r,dat,N,tmax)
print(real_ML)

LL(log(B),log(Phi),log(Y_0),0.23,dat,N,tmax)

LLhist <- numeric(0)
fit1 <- mle2(LL
        , start = list(log_B=log(B)
                     , log_Phi=log(Phi)
                     , logY_0=log(Y_0)
                     , r=r)
        , data = list(dat=dat
                    , N=N
                    , tmax=tmax
                    , debug = F
                    , debug_plot = F)
        , control = list(maxit=10000
                         ### parscale??
                       #, parscale = c(log(B), log(Phi), log(Y_0), r)
                         )
        , method = "Nelder-Mead"
        , skip.hessian = TRUE  ## TRUE to skip Hessian calculation ...
          )

## ?? not getting Re(ev) error any more?
print(real_ML)
print(-1*logLik(fit1))

coef(fit1)
true_param

### Testing with parameters away from real value

### Disturb B
# param <- list(log_B=log(0.01), log_Phi=log(Phi), logY_0=log(Y_0), r=r)
# param <- list(log_B=log(0.1), log_Phi=log(Phi), logY_0=log(Y_0), r=r)

## Identify true_param pretty well after shift to logistic.
## Hessian still not work
## Converge problem does not repeat for t=59
## When t=39 fit1 convergence failure: code=10 (degenerate Nelder-Mead simplex)
## However, fit2 works without error for t=39: conv problem is not necessary caused by initial values

### Disturb Phi
# param <- list(log_B=log(B), log_Phi=log(Phi+20), logY_0=log(Y_0), r=r)
# param <- list(log_B=log(B), log_Phi=log(Phi-15), logY_0=log(Y_0), r=r)
## Identify true_param pretty well after shift to logistic.
## Converge problem does not repeat for t=59, t=39

### Disturb Y_0
# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0+2e-4), r=r)
## Identify true_param pretty well after shift to logistic.
## Converge problem does not repeat for t=59, t=39

# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0-5e-5), r=r)
# function cannot be evaluated at initial parameters

# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0-3e-5), r=r)
# Works now

### Disturb r
# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r+0.2)
## function cannot be evaluated at initial parameters

param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r+0.1)
## t=39 not converging well, close but not enough log_lik=-547 while real param log_lik=-469
## increase maxit does not help
## t=59 function cannot be evaluated at initial parameters

param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r+0.05)
## indentify the init_param well

param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r-0.2)
## function cannot be evaluated at initial parameters
param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r-0.02)
## identify the init_param well
## sensitive to r


fit2 <- do.call(mle2,list(LL
                  , start = param
                  , data = list(dat=dat
                                , N=N
                                , tmax=tmax
                                , debug = T
                                , debug_plot = F)
                  , control = list(maxit=15000
                                   ### parscale??
                                   #, parscale = c(log(B), log(Phi), log(Y_0), r)
                  )
                  , method = "Nelder-Mead"
                  , skip.hessian = TRUE  ## TRUE to skip Hessian calculation ...
))

print(real_ML)
print(-1*logLik(fit2))
#print(fit2)
#summary(fit2)
param
coef(fit2)
## init_param



## re-do Hessian calculation with optimHess() ...
fix_hessian <- function(fit) {
    ## construct vectorized log-likelihood function
    lfun <- function(p) {
        do.call(c(as.list(p), fit@data), what = fit@minuslogl)
    }
    hh <- optimHess(coef(fit), fn = lfun)
    fit@vcov <- solve(hh)
    return(fit)
}

fit2H <- fix_hessian(fit2)
summary(fit2H)

## now try optimHess to see why we get NA values ...


fit2@details$hessian
## hmm, we get a finite hessian from this ...

quit()
hh <- optimHess(coef(fit1), fn = lfun)
vv <- solve(hh)
print(cov2cor(vv))
print(sdvec <- sqrt(diag(vv)))

## mle2 uses numDeriv::hessian() internally instead of optimHess() ...
numDeriv::hessian(lfun, coef(fit1))

warnings()

## still not sure why it's so hard to get a valid Hessian ...

