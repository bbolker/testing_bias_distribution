# Sys.setenv(LANG = "en")

library(shellpipes)
rpcall("const.fixed.pois.fit.Rout fixed.pois.fit.R const.pois.data.rda")
suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(macpan2))
options(macpan2_verbose = FALSE)

library(tidyr)
#library(ggplot2); theme_set(theme_bw())
loadEnvironments()

### test mechanism with True Value
S <- dat_all[which(dat_all$time==tmin-1 & dat_all$matrix=="S"),]$value
I <- dat_all[which(dat_all$time==tmin-1 & dat_all$matrix=="I"),]$value

sp_list <-tibble::lst(  beta=beta
                      , gamma=gamma
                      , N
                      , h=h
                      , w0=w0
                      , S=S
                      , I=I)


### Change simulation sir model

### To calculate Blik and T_B in simulator, we need data OT, OB as input
### Treat OP, OT as time varying parameter in the simulator
### Follwoing https://canmod.github.io/macpan2/articles/time_varying_parameters.html
### ???? Is there a better way to do this in MacPan2?

pos_changepoints = c(0,dat_fit$time-tmin+1)
neg_changepoints = c(0,dat_fit$time-tmin+1)
pos_values = c(0,dat_fit$OPos)
neg_values = c(0,dat_fit$ONeg)

### sqrt is removed from engine_functions now, cran has the update
### but not applied to website documents https://canmod.github.io/macpan2/reference/engine_functions.html

### Fitting model
(mc_sir
  |> mp_tmb_update(
    default = sp_list
  )
  |> mp_tmb_insert_backtrans(variables = c("beta","gamma","h","S","I"), mp_log)
  |> mp_tmb_insert(
    phase = "during"
    , at = Inf
    , expressions = list(
        pY ~ I/N                          ## Prevalence based on SIR
      , Phi ~ exp(-h)
      , pos ~ time_var(pos_values, pos_changepoints)
      , neg ~ time_var(neg_values, neg_changepoints)
      #, B_lik ~ 1/(2*N*Phi)*(((N-neg)*Phi+N-pos) - (((N-neg)*Phi+N-pos)^2-4*N*Phi*(N-pos-neg))^(1/2))
      , B_lik ~ exp(-w0)
      , T_B ~ 1 - B_lik
      , T_Y ~ 1 - Phi*B_lik
      # , T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
      # , T_pos ~ pY*T_Y/T_prop             ## Expected test positivity
      # , OT ~ rbinom(N,T_prop)
      # , OP ~ rbinom(OT,T_pos)
      , ONeg ~ T_B*(N-I)
      , OPos ~ T_Y*I
    )
    , default = list(  pos_values = pos_values
                     , neg_values = neg_values
    )
    , integers = list(  pos_changepoints = pos_changepoints
                      , neg_changepoints = neg_changepoints
                      )
  )
  |> mp_tmb_update(
    phase = "before"
    , at = 6
    , expressions = list(
      R ~ N - S - I
    )
  )
  ## |> mp_tmb_delete(phase = "before", at = Inf, default = c("beta","gamma","I","T_B","T_Y"))
) -> sir_sim

fit_pars <- c("log_beta", "log_gamma", "log_h", "log_w0" ,"log_S","log_I")
calibrator <- mp_tmb_calibrator(
    sir_sim
  , data = dat
  , traj = c("ONeg", "OPos", "T_B")
  , par = fit_pars,
  , default = list(N = N
                 , R = 0

    )
)
## modify likelihood function (eventually we'll have mp_binom() so we can do this
## when defining the calibrator)
calibrator$simulator$replace$obj_fn(~ - sum(dpois(obs_OPos, sim_OPos)) 
                                      - sum(dpois(obs_ONeg, sim_ONeg))
                                    )

fit<-mp_optimize(calibrator)
names(fit$par)<-fit_pars
#print(fit)

fit_bk<-c(exp(fit$par[1]),exp(fit$par[2]),N=N,exp(fit$par[3]),exp(fit$par[4]),exp(fit$par[5]),exp(fit$par[6]))

names(fit_bk)<-names(sp_list)
fit_bklist<-as.list(append(fit_bk,fit$par))

print(fit_bk)
print(unlist(sp_list))

print(fit$objective)

print(fit)

obj<-mp_tmb(calibrator)

library(numDeriv)
obj$fn(fit$par)

x<-jacobian(obj$gr,fit$par)
eigen(x)

saveEnvironment()