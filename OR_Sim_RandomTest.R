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
set.seed(13519)

## Random parameter values

# number of tests
n=100

### Constants for all tests
N <- 1e6         ## pop size
tmax <- 39       ## max simulation time/number of observations

### uninfected testing prob, from 0 to .25
T_B <- runif(2*n,0,0.25)
# 2n: Generate both true values and initial values for fitting

param_mat <- (expand.grid(T_B=T_B)
              %>% as_tibble()
              # Odds ratio part
              ## infected testing prob, must larger than T_B, smaller than 1
              %>% mutate(T_Y = runif(2*n,T_B,1))
              %>% mutate(B= T_B/(1-T_B))
              %>% mutate(Phi = (T_Y/(1-T_Y))/B)
              %>% mutate(log_B = log(B))
              %>% mutate(log_Phi = log(Phi))
              
              # Infection dynamic part
              %>% mutate(Y_0 = runif(2*n,0,25e-3))
              %>% mutate(log_Y0 = log(Y_0))
              ## initial value from 0 to .025 of population
              
              # logistic growth rate
              %>% mutate(r = log(2)/runif(2*n,1,5))
              ## random doubling time from 1 to 5
)

## separate simulation and fit values
param_true <- dplyr::slice(param_mat,1:n) 

param_fit <- dplyr::slice(param_mat,(n+1):(2*n))


## Simulate the observed data
dat_func <- function(param_vec, tmax, N){
  T_B <- param_vec$T_B
  T_Y <- param_vec$T_Y
  r <- param_vec$r
  Y_0 <- param_vec$Y_0
  t <- c(0:tmax)
  dat <- tibble(t=t
                , pY = 1/(1+(1/Y_0-1)*exp(-r*t))         ## Prevalence based on Logistic growth
                , T_prop = (1-pY)*T_B+pY*T_Y             ## Expected test proportion
                , pos = pY*T_Y/T_prop                    ## Expected test positivity
                , OT = rbinom(t,N,T_prop)                ## Observed number of test
                , OP = rbinom(t,OT,pos)                  ## Observed number of positive test
  )
  return(dat)
}

### function to calculate negative log-likelihood:
LL <- function(log_B, log_Phi, logY_0, r, dat, tmax, N, debug = FALSE) {
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

### fitting procedure
fit_proc <- function(dat,param_fit,tmax,N,debug=F){
  log_B <-param_fit$log_B
  log_Phi <- param_fit$log_Phi
  logY_0 <- param_fit$log_Y0
  r <- param_fit$r

  param <- list(log_B=log_B, log_Phi=log_Phi, logY_0=logY_0, r=r)
  fit <- mle2(LL
              , start = param
              , data = list(dat=dat
                            , N=N
                            , tmax=tmax
                            , debug = debug)
              , control = list(maxit=15000, reltol = 1e-10)
              , method = "Nelder-Mead"
  )
  return(fit)
}

true_logLik <- c(0)
fit_logLik <- c(0)
fit_logB <- c(0)
fit_logPhi <- c(0)
fit_logY0 <- c(0)
fit_r <- c(0)
inCI <- c(0)

results_list <- list(0)

for (i in c(1:n)) {
  dat <- dat_func(param_true[i,], tmax, N)
  true_logLik[i] <- LL(param_true[i,]$log_B,param_true[i,]$log_Phi,param_true[i,]$log_Y0,param_true[i,]$r,dat,tmax,N)

  fit <- fit_proc(dat,param_fit[i,],tmax,N)
  fit_logLik[i] <- -logLik(fit)
  fit_logB[i] <- coef(fit)[1]
  fit_logPhi[i] <- coef(fit)[2]
  fit_logY0[i] <- coef(fit)[3]
  fit_r[i] <- coef(fit)[4]

  true_value <- c(param_true[i,]$log_B,param_true[i,]$log_Phi,param_true[i,]$log_Y0,param_true[i,]$r)
  
  results<- tidy(fit, conf.int = TRUE) |>
    full_join(data.frame(term = tidy(fit)$term, true.value = true_value),
              by = "term") |>
    select(term, estimate, true.value, conf.low, conf.high)
  results_list[[i]]<-results
  inCI[i] <- min(results$true.value<=results$conf.high & results$true.value>=results$conf.low)
}

inCI
length(which(inCI==1))
length(which(is.na(inCI)))

NA_case <- which(is.na(inCI))
Zero_case <- which(inCI==0)
NA_case
Zero_case

### Checking cases
case <- 90

dat <- dat_func(param_true[case,], tmax, N)
dat

# matplot(dat$t, dat[,-1], type = "l", log = "y")
# legend("center", col = 1:4, lty = 1:4,
#        legend = names(dat)[-1])
# fit<-fit_proc(dat, param_fit[case,], tmax, N,debug=T)

select(param_true[case,],-log_Y0)
select(param_fit[case,],-log_Y0)
# 
# -logLik(fit)
# true_loglik
# 
# param_true[case,]
fit_logLik[case]
true_logLik[case]

true_value <- c(param_true[case,]$log_B,param_true[case,]$log_Phi,param_true[case,]$log_Y0,param_true[case,]$r)

## one way to present results ...
results_case <- results_list[[case]]

# results_case
inCI_case <- min(results_case$true.value<=results_case$conf.high & results_case$true.value>=results_case$conf.low)
inCI_case
ggplot(results_case, aes(y = term)) +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(x=true.value), colour = "red") +
  facet_wrap(~term, ncol = 1, scale  = "free")
