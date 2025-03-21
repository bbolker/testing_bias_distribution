# Sys.setenv(LANG = "en")
# remotes::install_github("bbolker/bbmle")

## for now we need a patched version of macpan2
## remotes::install_github("canmod/macpan2", ref = "dbinom2")

### ??? p_simulator dependence of "DEoptim" 
# install.packages("DEoptim")

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(macpan2))
suppressPackageStartupMessages(library(bbmle))
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(broom)
library(broom.mixed)
# library(DEoptim)

options(macpan2_verbose = FALSE)

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
tmin <- 30
tmax <- 80      ## max simulation time (about first half of logis)
# tmax <- 59     ## max simulation time (end of logis)
t <- c(tmin:tmax)
pts <- length(t) ## number of time points

## BMB: use self-naming list from tibble pkg
true_param <- tibble::lst(log_B=log(B), log_Phi=log(Phi), logY_0=log(Y_0), beta, gamma)

tp_list <-tibble::lst(beta, gamma, N, T_Y, T_B
              , I = NY_0
              , R = 0
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
    , outputs = c("pY","T_prop","pos","OT","OP","I")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  |> filter(time>=tmin)
) -> dat

dat
# dat$time <- dat$time-tmin+1

(dat
  |> pivot_wider(names_from = matrix,values_from = value)
) |> print(n=pts)

print(ggplot(dat)
 	+ aes(time, value, color=matrix)
 	+ geom_line()
 	+ scale_y_log10()
)

### Calibrator in macpan
## initial values for simulation
sp_list <-tibble::lst(beta, gamma, N, T_Y=T_Y, T_B
            , I = NY_0
            , R = 0
)

sir_sim <- mp_tmb_update(sir,
  default = sp_list
)
sir_sim |> mp_expand()

fit_pars <- c("beta", "gamma", "I", "T_Y", "T_B")
calibrator <- mp_tmb_calibrator(
    sir_sim
  , data = dat
  , traj = c("OT", "OP", "T_prop", "pos", "pY", "I")
  , par = fit_pars,
  , default = list(N = N
                 , R = 0
    )
  , time = mp_sim_bounds(1,tmax,"steps")
)
## modify likelihood function (eventually we'll have mp_binom() so we can do this
## when defining the calibrator)
calibrator$simulator$replace$obj_fn(~ - sum(dbinom(obs_OT, N, sim_T_prop)) - sum(dbinom(obs_OP, obs_OT, sim_pos)))

print(calibrator)

# fit_tp <- mp_optimize(calibrator)

fit_tp <- mp_optimize(calibrator)

cmp_par <- function(fit, true = unlist(tp_list),
                    pars = fit_pars) {
    data.frame(est = setNames(fit$par, pars),
               true = true[pars])
}

cmp_par(fit_tp)

## BMB: what was this for??                    , control=list(iter.max=10000, eval.max=10000))

print(fit_tp)

fit_loglik <- fit_tp$objective
print(fit_loglik)

## optim fit?

# fit_tp_optim <- mp_optimize(calibrator, "optim" ,method ="Nelder-Mead",control=list(maxit=10000, reltol = 1e-10))
# fit_tp_optim$convergence
# fit_loglik_optim<-fit_tp_optim$value
# print(fit_loglik_optim)


#mp_tmb_coef(calibrator)

# tmin and tmax would be a factor for fitting? 
# Could this degree of freedom shift to I_0 and N?

# mp_sim_offset not available!
# https://github.com/canmod/macpan2/blob/main/R/mp_tmb_calibrator.R
# https://github.com/canmod/macpan2/blob/main/man/mp_sim_offset.Rd
# ??? mp_sim_bounds()
# mp_sim_bounds(tmin,tmax,"steps")
# ??? mp_cal_time()
# need to change documentation! 
# ??? time argument of mp_tmb_calibrator

fit_tp_traj <- mp_trajectory(calibrator)
# fit_tp_traj
ggplot(fit_tp_traj, aes(time, value, color=matrix)) +
    geom_line() +
    scale_y_log10() +
    geom_point(data = dat, aes(x = time),size=0.8)

fit_tp_result <- (mp_tmb_coef(calibrator,conf.int = TRUE) 
                  |> select(-c("term", "type","row","col","std.error"))
                  |> cbind(true_value=unlist(tp_list)[fit_pars])
                  )

print(fit_tp_result)

## one way to present results ...
# results <- tidy(fit2, conf.int = TRUE) |> full_join(data.frame(term = names(true_param), true.value = true_param),
#               by = "term") |>
#     select(term, estimate, true.value, conf.low, conf.high)


## one way to show the results ...
# knitr::kable(results, digits = 3)


## or graphically ...
## (results are too precise, and range among true values is too large,
## to be able to see the confidence intervals if we plot everything on
## the same scale, so divide into separately scaled facets)
ggplot(fit_tp_result, aes(y = mat)) +
    geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
    geom_point(aes(x=true_value), colour = "red") +
    facet_wrap(~mat, ncol = 1, scale  = "free")

### does not work with beta+0.2
## false convergence (8)

p0 <- unlist(tp_list)[fit_pars]
tmb_obj <- mp_tmb(calibrator)
## nlminb and optim have different argument names but args 1, 2, 3 are
## starting value, objective function, gradient fn in both cases
(pars0 <- with(tmb_obj, nlminb(p0, fn, gr))$par)
(pars1 <- with(tmb_obj, optim (p0, fn, gr, method = "BFGS"))$par)
(pars2 <- with(tmb_obj, optim (p0+c(0.2, 0, 0, 0, 0), fn, gr, method = "BFGS"))$par)
(pars3 <- with(tmb_obj, nlminb(p0+c(0.2, 0, 0, 0, 0), fn, gr))$par)

(pmat <- cbind(pars0, pars1, pars2, pars3))
pmat - pmat[,1]

