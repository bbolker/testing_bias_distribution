# Sys.setenv(LANG = "en")
# remotes::install_github("bbolker/bbmle")

library(dplyr)
library(macpan2)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(broom)
source("mle2_tidy.R")

# Set Seeds
set.seed(13519)

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

# r <- log(2)/3    ## growth rate (doubling time = 3)
beta <- 0.25
gamma <- 0.1
tmin <- 50
tmax <- 80      ## max simulation time (about first half of logis)
# tmax <- 59     ## max simulation time (end of logis)
t <- c(tmin:tmax)
pts <- length(t) ## number of time points

true_param <- c("log_B"=log(B),"log_Phi"=log(Phi),"logY_0"=log(Y_0),"beta"=beta, "gamma"=gamma)

tp_list <- list(beta = beta
              , gamma = gamma
              , N = N
              , I = NY_0
              , R = 0
              , T_Y = T_Y
              , T_B = T_B
              )

### SIR from macpan
mc_sir <- mp_tmb_library("starter_models","sir", package = "macpan2")

(mc_sir 
  |> mp_tmb_update(
    default = tp_list
    )
  |> mp_tmb_insert(
      phase = "during"
    , at = Inf
    , expressions = list(
          pY ~ I/N                          ## Prevalence based on SIR
        , T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
        , pos ~ pY*T_Y/T_prop               ## Expected test positivity
        , OT ~ rbinom(N,T_prop)
        , OP ~ rbinom(OT,pos)
        )
    )
  )-> sir

# sir |> mp_expand()

(sir
  |> mp_simulator(
      time_steps = tmax
    , outputs = c("pY","T_prop","pos","OT","OP")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  |> filter(time>=tmin)
) -> dat

(dat
  |> pivot_wider(names_from = matrix,values_from = value)
) |> print(n=pts)

print(ggplot(dat)
 	+ aes(time, value, color=matrix)
 	+ geom_line()
 	+ scale_y_log10()
)

### Calibrator in macpan
# initial values for simulation
sp_list<-list(beta = beta
            , gamma = gamma
            , N = N
            , I = NY_0
            , R = 0
            , T_Y = T_Y
            , T_B = T_B
            )

sir_sim <- mp_tmb_update(sir,
  default = sp_list
)

calibrator <- mp_tmb_calibrator(
    sir_sim
  , data = dat
  , traj = c("OT","OP")
  , par = c("beta","gamma","I","T_Y","T_B")
  , default = list(N = N
                 , R = 0
    )
  )

print(calibrator)
### No mp_bin, just mp_neg_bin



# ### function to calculate negative log-likelihood:
# LL <- function(log_B, log_Phi, logY_0, beta, gamma
#                , dat, N, tmin ,tmax
#                , debug = FALSE
#                #,debug_plot = FALSE, plot_sleep = 1
#                ) {
#   Y_0 <- exp(logY_0)
#   NY_0 <- N*Y_0
#   B <- exp(log_B)
#   Phi <- exp(log_Phi)
#   T_B <- B/(1+B)
#   T_Y <- B*Phi/(1+B*Phi)
#   
#   t <- c(tmin:tmax)
#   pts <- length(t)
#   
#   beta <- beta
#   gamma <- gamma
#   
#   param_list <- list(  beta = beta
#                        , gamma = gamma
#                        , N = N
#                        , I = NY_0
#                        , R = 0
#                        , T_Y = T_Y
#                        , T_B = T_B
#   )
#   
#   (sir
#     |> mp_tmb_update(default = param_list)
#     |> mp_simulator(
#        time_steps = tmax
#       ,outputs = c("pY","T_prop","pos")
#       ) 
#     |> mp_trajectory()
#     |> dplyr::select(-c(row, col)) 
#     |> pivot_wider(names_from = matrix,values_from = value)
#     |> mutate(OT = rbinom(tmax,N,T_prop))
#     |> mutate(OP = rbinom(tmax,OT,pos))
#     |> dplyr::slice(tmin:tmax)
#     ) -> sim
#   ObsTest_nll <- -sum(dbinom(dat$OT, N, sim$T_prop, log = TRUE))
#   ObsPos_nll <- -sum(dbinom(dat$OP, dat$OT, sim$pos, log = TRUE))
#   out <- ObsTest_nll + ObsPos_nll
#   if (debug) {
#     cat(B, Phi, NY_0, beta, gamma, ObsTest_nll, ObsPos_nll, out, "\n")
#   }
#   return(out)
# }
# 
# real_ML <- LL(log(B),log(Phi),log(Y_0),beta,gamma,dat,N,tmin,tmax)
# print(real_ML)
# 
# LL(log(B),log(Phi),log(Y_0)+0.05,0.25,0.10,dat,N,tmin,tmax)
# 
# fit1 <- mle2(LL
#              , start = list(log_B=log(B)
#                             , log_Phi=log(Phi)
#                             , logY_0=log(Y_0)
#                             , beta=beta
#                             , gamma=gamma
#                             )
#              , data = list(dat=dat
#                            , N=N
#                            , tmin=tmin
#                            , tmax=tmax
#                            , debug = FALSE
#                            , debug_plot = FALSE
#                            )
#         , control = list(maxit=10000
#                          # parscale??
#                          #, parscale = c(log(B), log(Phi), log(Y_0), r)
#                          )
#         , method = "Nelder-Mead"
#         , hessian.method = "optimHess"
#         , skip.hessian = FALSE  ## TRUE to skip Hessian calculation ...
#         )
# 
# print(real_ML)
# print(-1*logLik(fit1))
# 
# coef(fit1)
# true_param
# 
# fit1@details$hessian
# ### This robust method provide an finite Hessian!

### Disturb B
# param <- list(log_B=log(0.01), log_Phi=log(Phi), logY_0=log(Y_0), beta=beta, gamma=gamma)
# param <- list(log_B=log(0.2), log_Phi=log(Phi), logY_0=log(Y_0), beta=beta, gamma=gamma)

## Identify init_param pretty well after shift to logistic
## Allowing wider parameter space
## Hessian works now

### Disturb Phi
param <- list(log_B=log(B), log_Phi=log(Phi+50), logY_0=log(Y_0), beta=beta, gamma=gamma)
# param <- list(log_B=log(B), log_Phi=log(Phi-20), logY_0=log(Y_0), beta=beta, gamma=gamma)

## Identify init_param pretty well after shift to logistic.
## Hessian works now, takes some time

### Disturb Y_0
# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0+2e-4), beta=beta, gamma=gamma)
## Identify init_param pretty well after shift to logistic.
## Converge problem does not repeat for t=59, t=39

# param <- list(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0-5e-5), beta=beta, gamma=gamma)
## Works for smaller Y_0 value now

### Disturb beta

### Disturb gamma

## profiling showed that we can get a slightly better fit ...
## decreasing tolerance avoids that problem
fit2 <- mle2(LL
           , start = param
           , data = list(dat=dat
                       , N=N
                       , tmin=tmin
                       , tmax=tmax
                       , debug = TRUE)
           , control = list(  maxit=10000
                            , reltol = 1e-10
                            )
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

# NY_0_fit2 <- N*exp(coef(fit2)[3])
# beta_fit2 <- coef(fit2)[4]
# gamma_fit2 <- coef(fit2)[5]


## one way to present results ...
results <- tidy(fit2, conf.int = TRUE) |> full_join(data.frame(term = names(true_param), true.value = true_param),
              by = "term") |>
    select(term, estimate, true.value, conf.low, conf.high)

# # ## a little slow (6 seconds)
# system.time(
#     results_prof <- tidy(fit2, conf.int = TRUE, conf.method = "spline")
# )
# 
# ## very little difference in this case (although CIs are narrow anyway)
# results_prof$conf.low-results$conf.low
# results_prof$conf.high-results$conf.high

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

