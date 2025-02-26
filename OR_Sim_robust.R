# Sys.setenv(LANG = "en")
# remotes::install_github("bbolker/bbmle")

library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(broom)
source("mle2_tidy.R")

# Set Seeds
# set.seed(13519)

## Initial true values:
T_B <- 0.04               ## uninfected testing prob
T_Y <- 0.5                ## infected testing prob
B <- T_B/(1-T_B)          ## baseline odds of testing
Phi <- (T_Y/(1-T_Y))/B    ## inf vs uninf testing odds ratio

print(B)
print(Phi)

Y_0 <- 1e-4      ## initial prevalence
N <- 1e6         ## pop size
NY_0 <- N*Y_0    ## initial number infected

r <- log(2)/3    ## growth rate (doubling time = 3)
tmax <- 39       ## max simulation time (about first half of logis)
# tmax <- 59     ## max simulation time (end of logis)
t <- c(0:tmax)
pts <- length(t) ## number of time points

true_param <- c("log_B"=log(B),"log_Phi"=log(Phi),"logY_0"=log(Y_0),"r"=r)


## Simulate the data
dat <- tibble(t=t
	## , pY = pmin(Y_0*exp(r*t), 1)          ## Exponential growth
	, pY = 1/(1+(1/Y_0-1)*exp(-r*t))         ## Prevalence based on Logistic growth
	, T_prop = (1-pY)*T_B+pY*T_Y             ## Expected test proportion
	, pos = pY*T_Y/T_prop                    ## Expected test positivity
	, OT = rbinom(t,N,T_prop)                ## Observed number of test
	, OP = rbinom(t,OT,pos)                  ## Observed number of positive test
)

# print(dat,n=60)
matplot(dat$t, dat[,c(-1,-3)], type = "l", log = "y")
legend("center", col = 1:4, lty = 1:4,
       legend = names(dat)[c(-1,-3)])

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
LL <- function(log_B, log_Phi, logY_0, r, dat, N, tmax, debug = FALSE,
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
                  ## , pY = pmin(Y_0*exp(r*t), 1)          ## Exponential growth
                  , pY = 1/(1+(1/Y_0-1)*exp(-r*t))         ## Prevalence based on Logistic growth
                  , T_prop = (1-pY)*T_B+pY*T_Y             ## Expected test proportion
                  , pos = pY*T_Y/T_prop                    ## Expected test positivity
    )
  # if(max(sim$pY) == 1 || any(sim$NY<dat$posTests) || any((N-sim$NY)<dat$negTests) || any(N<sim$NY)) return(NA)

  # if (any(sim$NY<dat$posTests)) {
  #     cat("Underestimated infected population, pos tests > infected population", "\n")
  # }
  # if (any((N-sim$NY)<dat$negTests)) {
  #     cat("Overestimated infected population, neg tests > uninfected population", "\n")
  # }
  ObsTest_nll <- -sum(dbinom(dat$OT, N, sim$T_prop, log = TRUE))
  ObsPos_nll <- -sum(dbinom(dat$OP, dat$OT, sim$pos, log = TRUE))
  out <- ObsTest_nll + ObsPos_nll
  if (debug) {
      cat(B, Phi, logY_0, r, ObsTest_nll, ObsPos_nll,
          out, "\n")
  }
  # if (debug_plot) {
  #     par(mfrow= c(1,2), las = 1)
  #     ylim <- range(c(dat$posTests, dat$negTests,
  #                            sim$NY*T_Y, (N-sim$NY)*T_B))
  #     matplot(dat$t, dat[c("posTests", "negTests")], type = "p",
  #             pch = 1:2, log = "y",
  #             ylim = ylim)
  #     matlines(dat$t, cbind(sim$NY*T_Y, (N-sim$NY)*T_B))
  #     LLhist <<- c(LLhist, out)
  #     plot(LLhist - min(LLhist) + 1e-3, type = "b", log = "y")
  #     Sys.sleep(plot_sleep)
  # }
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
        , hessian.method = "optimHess"
        , skip.hessian = FALSE  ## TRUE to skip Hessian calculation ...
          )

print(real_ML)
print(-1*logLik(fit1))

coef(fit1)
true_param

fit1@details$hessian
### This robust method provide an finite Hessian!

### Disturb B
# param <- list(log_B=log(0.01), log_Phi=log(Phi), logY_0=log(Y_0), r=r)
# param <- list(log_B=log(0.2), log_Phi=log(Phi), logY_0=log(Y_0), r=r)

## Identify init_param pretty well after shift to logistic
## Allowing wider parameter space
## Hessian works now

### Disturb Phi
# param <- list(log_B=log(B), log_Phi=log(Phi+50), logY_0=log(Y_0), r=r)
# param <- list(log_B=log(B), log_Phi=log(Phi-20), logY_0=log(Y_0), r=r)

## Identify init_param pretty well after shift to logistic.
## Hessian works now, takes some time

### Disturb Y_0
# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0+2e-4), r=r)
## Identify init_param pretty well after shift to logistic.
## Converge problem does not repeat for t=59, t=39

# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0-5e-5), r=r)
## Works for smaller Y_0 value now

### Disturb r
# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r+0.2)
## Works for larger r value now

param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), r=r-0.2)
## Works for lower r value now
## sensitive to r

## profiling showed that we can get a slightly better fit ...
## decreasing tolerance avoids that problem
fit2 <- mle2(LL
           , start = param
           , data = list(dat=dat
                       , N=N
                       , tmax=tmax)
           , control = list(maxit=15000, reltol = 1e-10)
           , method = "Nelder-Mead"
)

print(real_ML)
print(-1*logLik(fit2))
print(-1*logLik(fit1))
# print(fit2)

#param
coef(fit2)
true_param

summary(fit2)
vcov(fit2)

## one way to present results ...

results <- tidy(fit2, conf.int = TRUE) |>
    full_join(data.frame(term = names(true_param), true.value = true_param),
              by = "term") |>
    select(term, estimate, true.value, conf.low, conf.high)

# ## a little slow (6 seconds)
# system.time(
#     results_prof <- tidy(fit2, conf.int = TRUE, conf.method = "spline")
# )

## very little difference in this case (although CIs are narrow anyway)
results_prof$conf.low-results$conf.low
results_prof$conf.high-results$conf.high

## one way to show the results ...
knitr::kable(results, digits = 3)

## or graphically ...

## (results are too precise, and range among true values is too large,
##  to be able to see the confidence intervals if we plot everything on
## the same scale, so divide into separately scaled facets)
ggplot(results, aes(y = term)) +
    geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
    geom_point(aes(x=true.value), colour = "red") +
    facet_wrap(~term, ncol = 1, scale  = "free")

## we would like to compute profile confidence intervals, but this is slightly
## problematic
# pp0 <- profile(fit2)
# logLik(pp0)
# logLik(fit2)
# 
# cbind(coef(pp0), coef(fit2))


### Randomize initial parameter for fitting
TB_random <- runif(1,0,0.25)
TY_random <- runif(1,TB_random,1)
B_random <- TB_random/(1-TB_random)
logB_random <- log(B_random)
Phi_random <- (TY_random/(1-TY_random))/B
logPhi_random <- log(Phi_random)

Y0_random <- round(runif(1,0,5e-4),6)
logY0_random <- log(Y0_random)
r_random <- log(2)/runif(1,0,5)

param_rd_vec <- c("log_B"=logB_random,"log_Phi"=logPhi_random,"logY_0"=logY0_random,"r"=r_random)

param_rd <- list(log_B=logB_random, log_Phi=logPhi_random, logY_0=logY0_random, r=r_random)

fit3 <- mle2(LL
             , start = param_rd
             , data = list(dat=dat
                           , N=N
                           , tmax=tmax
                           , debug = T)
             , control = list(maxit=25000, reltol = 1e-10)
             , method = "Nelder-Mead"
)

print(real_ML)
print(-1*logLik(fit3))
vcov(fit3)
#param
coef(fit3)
true_param
param_rd_vec

## one way to present results ...
results3 <- tidy(fit3, conf.int = TRUE) |>
  full_join(data.frame(term = names(true_param), true.value = true_param),
            by = "term") |>
  select(term, estimate, true.value, conf.low, conf.high)
results3
## one way to show the results ...
# knitr::kable(results3, digits = 3)

## or graphically ...

## (results are too precise, and range among true values is too large,
##  to be able to see the confidence intervals if we plot everything on
## the same scale, so divide into separately scaled facets)
ggplot(results3, aes(y = term)) +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(x=true.value), colour = "red") +
  facet_wrap(~term, ncol = 1, scale  = "free")
