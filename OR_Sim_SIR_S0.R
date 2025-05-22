# Sys.setenv(LANG = "en")
## remotes::install_github("bbolker/bbmle")

## for now we need a patched version of macpan2
## remotes::install_github("canmod/macpan2", ref = "dbinom2")
## remotes::install_github("canmod/macpan2")

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
library(DEoptim)
library(nloptr)
## https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-gradient-based-optimization

options(macpan2_verbose = FALSE)

source("mle2_tidy.R")
source("spec_trans_par.R")

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
tmax <- 80      ## max simulation time 
# tmax <- 59     ## max simulation time 
t <- c(tmin:tmax)
pts <- length(t) ## number of time points

logit_trans <- function(x){
  log(x)-log(1-x)
}
logit_backtrans <- function(x){
  (1/(1 + exp(-x)))
}

## BMB: use self-naming list from tibble pkg
true_param <- tibble::lst(  log_B=log(B)
                          , log_Phi=log(Phi)
                          , logY_0=log(Y_0)
                          , logit_T_B=logit_trans(T_B)
                          , logit_T_Y=logit_trans(T_Y)
                          )

tp_list <-tibble::lst(beta, gamma, N, T_B, T_Y
              , I = NY_0
              , R = 0
)

### SIR from macpan
mc_sir <- mp_tmb_library("starter_models","sir", package = "macpan2")

(mc_sir
  |> mp_tmb_update(
    default = tp_list
    )
  |> mp_tmb_insert_backtrans(variables = c("beta","gamma","I"), mp_log)
  |> mp_tmb_insert_backtrans(variables = c("T_B","T_Y"), mp_logit)
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
  ## |> mp_tmb_delete(phase = "before", at = Inf, default = c("beta","gamma","I","T_B","T_Y"))
)->sir

sir |> mp_expand()

# sir |> mp_default()

(sir
  |> mp_simulator(
      time_steps = tmax
    , outputs = c("pY","T_prop","pos","OT","OP","I","S")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  
) -> dat_all

dat<- (dat_all|> filter(time>=tmin)
)

# dat
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

### Assume we know St
S <- dat_all[which(dat_all$time==tmin-1 & dat_all$matrix=="S"),]$value
hat_S <- S*(1-0.4)

### approximation of hat{I} inferred from first data point:
OT <- dat[which(dat$time==tmin & dat$matrix=="OT"),]$value
OP <- dat[which(dat$time==tmin & dat$matrix=="OP"),]$value

print(T_Y)
hat_T_Y <- T_Y+0.2

hat_T <- OT/N
hat_p <- OP/OT
hat_Y <- (hat_p*hat_T)/hat_T_Y

hat_I <- hat_Y*N

### In fitting, should check if S+I >=1
if(hat_S+hat_I>N){
  print("initial S+I value larger than N")
} else {"check"}


### What if we don't know S???
# it will only affect the scaling 
### How can we detect S???

sp_list <-tibble::lst(beta=beta+0.4, gamma=gamma+0.1, N, T_B, T_Y=hat_T_Y, S=hat_S, I=hat_I)

### Change simulation sir model
(mc_sir
  |> mp_tmb_update(
    default = sp_list
  )
  |> mp_tmb_insert_backtrans(variables = c("beta","gamma","S","I"), mp_log)
  |> mp_tmb_insert_backtrans(variables = c("T_B","T_Y"), mp_logit)
  |> mp_tmb_insert(
    phase = "during"
    , at = Inf
    , expressions = list(
      pY ~ I/N                            ## Prevalence based on SIR
      , T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
      , pos ~ pY*T_Y/T_prop               ## Expected test positivity
      , OT ~ rbinom(N,T_prop)
      , OP ~ rbinom(OT,pos)
    )
  )
  |> mp_tmb_update(
    phase = "before"
    , at = 7
    , expressions = list(
      R ~ N - S - I
    )
  )
  ## |> mp_tmb_delete(phase = "before", at = Inf, default = c("beta","gamma","I","T_B","T_Y"))
) -> sir_sim

# print(sir_sim) 

fit_pars <- c("log_beta", "log_gamma", "log_S","log_I", "logit_T_B", "logit_T_Y")
calibrator <- mp_tmb_calibrator(
    sir_sim
  , data = dat
  , traj = c("OT", "OP", "T_prop", "pos")
  , par = fit_pars,
  , default = list(N = N
                 , R = 0
    )
)
## modify likelihood function (eventually we'll have mp_binom() so we can do this
## when defining the calibrator)
calibrator$simulator$replace$obj_fn(~ - sum(dbinom(obs_OT, N, sim_T_prop)) - sum(dbinom(obs_OP, obs_OT, sim_pos)))

#calibrator|>print()

#mp_optimize(calibrator,optimizer = "optim", method = "BFGS")
#mp_optimize(calibrator)


fit<-mp_optimize(calibrator)
names(fit$par)<-fit_pars
print(fit)


fit_bk<-c(exp(fit$par[1]),exp(fit$par[2]),N=N,logit_backtrans(fit$par[5]),logit_backtrans(fit$par[6]),exp(fit$par[3]),exp(fit$par[4]))

names(fit_bk)<-names(sp_list)
fit_bklist<-as.list(append(fit_bk,fit$par))

# print(fit_bk)

## Look into initial values
test_list <- sp_list


(sir_sim
  |> mp_tmb_update(phase = "during", default = test_list)
  |> mp_tmb_insert_backtrans(variables = c("beta","gamma","S","I"), mp_log)
  |> mp_tmb_insert_backtrans(variables = c("T_B","T_Y"), mp_logit)
  |> mp_simulator(
    time_steps = tmax-tmin+1
    , outputs = c(  "OP"
                  , "OT"
                  , "pY"
                  )
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  # |> filter(time>=tmin)
) -> dat_sim


# (sir_sim
#   |> mp_tmb_update(phase = "during", default = fit_bklist)
# )

(sir_sim|> mp_tmb_update(phase = "during", default = fit_bklist)
        |> mp_simulator(
            time_steps = tmax-tmin+1
            , outputs = c(  "OP"
                          , "OT"
                          , "pY"
                          )
             ) 
)->sir_optim
sim_vals <-(sir_optim|> mp_trajectory()
              |> dplyr::select(-c(row, col))
             # |> filter(time>=tmin)
)
mp_default(sir_optim)


sir_optim
fit_bklist
sim_vals[1,]

# (dat_sim
#   |> pivot_wider(names_from = matrix,values_from = value)
# ) |> print(n=pts)
dat_tp<-filter(dat,matrix=="OP"|matrix=="OT"|matrix=="pY")

dat_sim$time<-dat_sim$time+tmin-1
sim_vals$time<-sim_vals$time+tmin-1

dat_sim <- cbind(dat_sim,model=rep("sim_init",length(dat_sim[,1])))
dat_tp <- cbind(dat_tp,model=rep("real",length(dat_tp[,1])))
sim_vals <- cbind(sim_vals,model=rep("optim",length(sim_vals[,1])))

dat_compare<-rbind(dat_sim,dat_tp,sim_vals)

fit_curve <- (ggplot(dat_compare)
      + aes(x=time, y=value, color=matrix, linetype = model)
      + geom_line()
      + geom_point(data=dat_tp,aes(shape="real"))
      + scale_y_log10()
      #+ scale_colour_manual(values = c("blue", "red", "black"))
      + scale_linetype_manual(values = c(1,2,3))
)
print(fit_curve)
ggsave("new_mech_curve.png",plot=fit_curve, path = "./pix", width=3200,height=1800,units="px")

### Obs: difference between beta+0.25 and beta+0.30 is if the tipping point is 
### contained in the initial simulation

### Obs: lower the tmin to include the tipping/peak point will fix the 
### indentifiability issue

### Conjecture: including the peak or have correct tendency make a difference in fitting 
### similar things happens when tmin=20 to beta+0.45 and beta+0.50
### A fixable issue??: make sure the initial value for testing at least have the same tendency with data




# ## NOT IDEMPOTENT, mutability issues,  etc .... ???
# par_trans <- c(beta = "log", gamma = "log", T_Y = "logit", T_B = "logit", I = "log")
# cal_trans_spec <- mp_trans_pars(sir_sim, par_trans)
# 
# calibrator_trans <- mp_tmb_calibrator(
#     cal_trans_spec
#   , data = dat
#   , traj = c("OT", "OP", "T_prop", "pos", "pY", "I")
#     ## use transformed names
#   , par = mk_par_names(par_trans),
#   , default = list(N = N
#                  , R = 0
#     )
#   , time = mp_sim_bounds(1,tmax,"steps")
# )
# calibrator_trans$simulator$replace$obj_fn(~ - sum(dbinom(obs_OT, N, sim_T_prop)) - sum(dbinom(obs_OP, obs_OT, sim_pos)))
# 
# ## modify?
# ## calibrator_trans <- mp_trans_args(calibrator, par_trans)


#' @param p parameter vector
#' @param off_par names or indices of parameters to modify
#' @param off_val values to offset specified parameters
#' @param cal macpan calibrator object
#' @param optimizer name (or symbol) of optimizer function
#' @param ret_val "all" to return full list, character vector to return subset of list
#' @param ... additional arguments (e.g. method for optim)
# my_opt <- function(p,
#                    off_par = NULL,
#                    off_val = 0,
#                    cal = calibrator,
#                    optimizer = "nlminb",
#                    ret_val = "par", ...) {
#     obj <- mp_tmb(cal)
#     if (!is.null(off_par)) p[off_par] <- p[off_par]+off_val
#     if (is.character(optimizer)) optimizer <- get(optimizer)
#     log_p <- log(p)
#     names(log_p) <- fit_pars
#     fit <- optimizer(log_p, obj$fn, obj$gr, ...)
#     if (ret_val == "all") return(fit)
#     return(fit[ret_val])
# }
# 
# 
# p0 <- unlist(tp_list)[c("beta","gamma","I","T_B","T_Y")]
# 
# my_opt(p0)
# 
# my_opt(p0, off_par = "beta", off_val = 0.3, ret_val = "objective")
# 
# offvec <- seq(-0.24, 0.6, by = 0.01)
# nllvec <- sapply(offvec, \(x) my_opt(p0, off_par = "beta", off_val = x, ret_val = "objective")[[1]])
# ## max value that's OK
# offvec[which(nllvec>1000)[1] -1 ]
# nllvec
# 
# nllvec2 <- sapply(offvec, \(x) my_opt(p0, off_par = "beta", off_val = x,
#                                       optimizer = "optim", method = "BFGS",
#                                       ret_val = "value")[[1]])
# nllvec2
# 
# 
# nlfun <- function(par, fn, gr, algorithm = "NLOPT_LD_LBFGS") {
#     ## hack around nloptr limitations
#     fn0 <- function(x) fn(x)
#     gr0 <- function(x) gr(x)
#     fit <- suppressWarnings(nloptr(x0 = par, eval_f = fn0, eval_grad_f = gr0, opts = list(algorithm = algorithm)))
#     return(list(par = fit$solution, value = fit$objective, convergence = fit$status))
# }
# 
# my_opt(p0, optimizer = "nlfun", ret_val = "all")
# 
# nllvec3 <- sapply(offvec, \(x) my_opt(p0, off_par = "beta", off_val = x,
#                                       optimizer = "nlfun",
#                                       ret_val = "value")[[1]])
# 
# # fit_tp <- mp_optimize(calibrator)
# 
# fit_tp <- mp_optimize(calibrator)
# 
# cmp_par <- function(fit, true = unlist(tp_list),
#                     pars = fit_pars) {
#     data.frame(est = setNames(fit$par, pars),
#                true = true[pars])
# }
# 
# cmp_par(fit_tp)

## BMB: what was this for??                    , control=list(iter.max=10000, eval.max=10000))

# print(fit_tp)
# 
# fit_loglik <- fit_tp$objective
# print(fit_loglik)

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
# ggplot(fit_tp_traj, aes(time, value, color=matrix)) +
#     geom_line() +
#     scale_y_log10() +
#     geom_point(data = dat, aes(x = time),size=0.8)
mp_tmb_coef(calibrator) 
dat_tp[1,]$value*N
sim_vals[1,]$value*N
fit_bklist$I

fit_tp_result <- (mp_tmb_coef(calibrator,conf.int = TRUE) 
                  |> select(-c("term", "type","row","col","std.error"))
                  |> cbind(true_value=c(beta=beta,gamma=gamma,S=S,I=dat_all[197,]$value,T_B=T_B,T_Y=T_Y)
)
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
fit_compare <- ggplot(fit_tp_result, aes(y = mat)) +
    geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
    geom_point(aes(x=true_value), colour = "red") +
    facet_wrap(~mat, ncol = 1, scale  = "free")
print(fit_compare)
ggsave("new_mech_fit.png",plot=fit_compare, path = "./pix", width=1800,height=3200,units="px")

### does not work with beta+0.2
## false convergence (8)

# p0 <- unlist(tp_list)[fit_pars]
# tmb_obj <- mp_tmb(calibrator)
# ## nlminb and optim have different argument names but args 1, 2, 3 are
# ## starting value, objective function, gradient fn in both cases
# (pars0 <- with(tmb_obj, nlminb(p0, fn, gr))$par)
# (pars1 <- with(tmb_obj, optim (p0, fn, gr, method = "BFGS"))$par)
# (pars2 <- with(tmb_obj, optim (p0+c(0.3, 0, 0, 0, 0), fn, gr, method = "BFGS"))$par)
# (pars3 <- with(tmb_obj, nlminb(p0+c(0.3, 0, 0, 0, 0), fn, gr))$par)
# 
# (pmat <- cbind(pars0, pars1, pars2, pars3))
# pmat - pmat[,1]
# 
