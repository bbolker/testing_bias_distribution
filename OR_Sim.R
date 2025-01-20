
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
tmax <- 19       ## max simulation time
t <- c(0:tmax)
pts <- length(t) ## number of time points

## Simulate the data
dat <- tibble(t=t
	, pY = pmin(Y_0*exp(r*t), 1)
	## , NY = rbinom(pts, N, pY)
	, NY = round(N*pY)
	, posTests = rbinom(pts, NY, T_Y)
	, negTests = rbinom(pts, N-NY, T_B)
)

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
                , pY = pmin(Y_0*exp(r*t), 1)
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

LL(log(0.05),log(Phi),log(Y_0),0.2,dat,N,tmax)

LLhist <- numeric(0)
fit1 <- mle2(LL
        , start = list(log_B=log(B)
                     , log_Phi=log(Phi)
                     , logY_0=log(Y_0)
                     , r=r)
        , data = list(dat=dat
                    , N=N
                    , tmax=tmax
                    , debug = FALSE)
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


### Testing with parameters away from real value
# param <- list(log_B=log(0.02), log_Phi=log(Phi), logY_0=log(Y_0), r=r)

param <- list(log_B=log(0.07), log_Phi=log(Phi), logY_0=log(Y_0), r=r)
## convergence failure: code=10 (degenerate Nelder-Mead simplex)

# param <- list(log_B=log(B), log_Phi=log(Phi+10), logY_0=log(Y_0), r=r)
## convergence failure: code=10 (degenerate Nelder-Mead simplex)

# param <- list(log_B=log(B), log_Phi=log(Phi-15), logY_0=log(Y_0), r=r)

# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0+1e-4), r=r)
## sensitive to Y_0: Phi is larger

fit2 <- do.call(mle2,list(LL
                  , start = param
                  , data = list(dat=dat
                                , N=N
                                , tmax=tmax
                                , debug = FALSE
                                , debug_plot = FALSE)
                  , control = list(maxit=10000
                                   ### parscale??
                                   #, parscale = c(log(B), log(Phi), log(Y_0), r)
                  )
                  , method = "Nelder-Mead"
                  , skip.hessian = TRUE  ## TRUE to skip Hessian calculation ...
))

print(real_ML)
print(-1*logLik(fit2))
print(fit2)
summary(fit2)


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


fit1@details$hessian
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

