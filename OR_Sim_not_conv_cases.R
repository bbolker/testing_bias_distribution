# Sys.setenv(LANG = "en")
# remotes::install_github("bbolker/bbmle")

library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(broom)
# source("mle2_tidy.R")

# Set Seeds
set.seed(13519)

## Random parameter values

## BMB: good to be consistent about <- or = for assignment
## see e.g. https://lintr.r-lib.org/
##  https://www.tidyverse.org/blog/2025/02/air/

# number of tests
n <- 1000

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
              ## BMB: you can chain mutate expressions
              %>% mutate(T_Y = runif(2*n,T_B,1),
                         B= T_B/(1-T_B),
                         Phi = (T_Y/(1-T_Y))/B,
                         log_B = log(B),
                         log_Phi = log(Phi),
                         
                         ## Infection dynamic part
                         Y_0 = runif(2*n,0,25e-3),
                         log_Y0 = log(Y_0),
                         ## initial value from 0 to .025 of population
              
                         ## logistic growth rate
                         r = log(2)/runif(2*n,1,5)
                         ## random doubling time from 1 to 5
                         )
)

## BMB: I probably wouldn't bother with dplyr::slice here
##  (just param_mat[1:n])
## separate simulation and fit values
param_true <- dplyr::slice(param_mat,1:n) 

## initial values for fitting
param_fit <- dplyr::slice(param_mat,(n+1):(2*n))

## Simulate the observed data
dat_func <- function(param_vec, tmax, N) {
  ##  BMB: could use with(param_vec, { ... } ) (at slight loss of debugging capability)
  T_B <- param_vec$T_B
  T_Y <- param_vec$T_Y
  r <- param_vec$r
  Y_0 <- param_vec$Y_0
  t <- c(0:tmax)
  dat <- tibble(t=t
              , pY = 1/(1+(1/Y_0-1)*exp(-r*t))  ## Prevalence based on Logistic growth
              , T_prop = (1-pY)*T_B+pY*T_Y      ## Expected test proportion
              , pos = pY*T_Y/T_prop             ## Expected test positivity
              , OT = rbinom(t,N,T_prop)         ## Observed number of tests
              , OP = rbinom(t,OT,pos)           ## Observed number of positive tests
  )
  return(dat)
}

## BMB: changed logY_0 to logY_0 for consistency (could also be
## log_Y0 or log_Y_0, but aim for consistency in any case)

### function to calculate negative log-likelihood:
LL <- function(log_B, log_Phi, log_Y0, r, dat, tmax, N, debug = FALSE,
               debug_plot = FALSE, plot_sleep = 1) {
    Y_0 <- exp(log_Y0)
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
  ObsTest_nll <- -sum(dbinom(dat$OT, N, sim$T_prop, log = TRUE))
  ObsPos_nll <- -sum(dbinom(dat$OP, dat$OT, sim$pos, log = TRUE))
  out <- ObsTest_nll + ObsPos_nll
  if (debug) {
      cat(B, Phi, log_Y0, r, ObsTest_nll, ObsPos_nll,
          out, "\n")
  }
  return(out)
}

### fitting procedure
fit_proc <- function(dat,param_fit,tmax,N,debug=FALSE){

  log_B <-param_fit$log_B
  log_Phi <- param_fit$log_Phi
  log_Y0 <- param_fit$log_Y0
  r <- param_fit$r

  ## BMB: why do you unpack param_fit and then re-pack it? Seems unnecessary
    
  param <- list(log_B=log_B, log_Phi=log_Phi, log_Y0=log_Y0, r=r)
    
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

## All 10 cases that fail to converge to true values
NA_case <- c(227,313,343,436,660,711,725,760,882,956)

## plot joint distributions of starting values and true values
plong <- function(x, lab = "true_val") {
    (x
        |> select(-c(log_Phi, log_B, log_Y0))
        |> mutate(run = seq(n()))
        |> pivot_longer(-run,
                        names_to = "param",
                        values_to = lab)
    )
}
comb <- full_join(plong(param_fit, "fit_val"),
          plong(param_true, "true_val"),
          by = c("run", "param"))

ggplot(comb, aes(fit_val, true_val)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~param, scale = "free") +
    geom_point(data = comb[NA_case,], col = "red", size = 3)


param_fit[NA_case,]
### Checking cases
case <- NA_case[1]

select(param_true[case,],-log_Y0)
select(param_fit[case,],-log_Y0)

## fit_logLik[case]
# true_logLik[case]

### fit again with debug on
dat <- dat_func(param_true[case,], tmax, N)
# fitting result might change due to randomness in data generating

# matplot(dat$t, dat[,-1], type = "l", log = "y")
# legend("center", col = 1:4, lty = 1:4,
#        legend = names(dat)[-1])

vars <- c("log_B", "log_Phi", "log_Y0","r")
true_logLik <- do.call(LL, c(param_true[case, vars], list(dat, tmax, N)))

fit_case<- fit_proc(dat, param_fit[case,], tmax, N,debug=TRUE)

fit_logLik <- -logLik(fit_case)

true_logLik
fit_logLik

true_value <- param_true[case, vars]

### error starts here due to ill-behaved or missing hessian
tt <- tidy(fit_case, conf.int = TRUE)
results_case <- tt |>
    full_join(data.frame(term = names(true_value),
                         true.value = unlist(true_value)),
              by = "term") |>
    select(term, estimate, true.value, conf.low, conf.high)

# results_case
# inCI_case <- min(results_case$true.value<=results_case$conf.high & results_case$true.value>=results_case$conf.low)
# inCI_case
# ggplot(results_case, aes(y = term)) +
#   geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
#   geom_point(aes(x=true.value), colour = "red") +
#   facet_wrap(~term, ncol = 1, scale  = "free")




