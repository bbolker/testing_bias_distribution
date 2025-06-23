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
# library(broom)
# library(broom.mixed)
# library(DEoptim)
# library(nloptr)
## https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-gradient-based-optimization

options(macpan2_verbose = FALSE)

my_sir_dir = file.path(getwd(), "sirs_seasonal")
# mp_model_starter("sir", my_sir_dir)
sir_season = mp_tmb_library(my_sir_dir)

(sir_season
  |> mp_tmb_update(  phase="before",
                   , default=list(  gamma=0.05
                                  , beta_high = 0.2
                                  , beta_low = 0.05
                                  , eta = 0.01 ) 
                   )
  |> mp_simulator(
      time_steps = 3000
    , outputs = c ("I","beta","S")
  )
  |> mp_trajectory()
  |> mutate(quantity = case_match(matrix
                                  , "I" ~ "Prevalence"
                                  , "beta" ~ "beta"
                                  , "S" ~ "S"
                                  )
            )
) -> naive_traj

(naive_traj|> ggplot() 
  + geom_line(aes(time, value)) 
  + facet_wrap(~ quantity, scales = "free")
  + theme_bw()
)

eq_SI<-naive_traj[which(naive_traj$time==365*5+1),]$value
I_0 <- eq_SI[1]
S_0 <- eq_SI[2]
R_0 <- 10000-S_0-I_0

(sir_season
  |> mp_tmb_update(  phase="before",
                   , default=list(  gamma=0.05
                                  , beta_high = 0.2
                                  , beta_low = 0.05
                                  , eta = 0.01
                                  , I = I_0
                                  , R = R_0
                                  ) 
                   )
  |> mp_simulator(
                    time_steps = 1000
                  , outputs = c ("I","beta","S")
                  )
  |> mp_trajectory()
  |> mutate(quantity = case_match(matrix
                                , "I" ~ "Prevalence"
                                , "beta" ~ "beta"
                                , "S" ~ "S"
                                )
            )
  |> ggplot() 
  + geom_line(aes(time, value)) 
  + facet_wrap(~ quantity, scales = "free")
  + theme_bw()
)
######################################


### Seasonal SIRS model constructed
### Simulate "ture" data

### Ture parameters
## Initial true values:
T_B <- 0.02               ## uninfected testing prob
T_Y <- 0.5                ## infected testing prob
B <- T_B/(1-T_B)          ## baseline odds of testing
Phi <- (T_Y/(1-T_Y))/B    ## inf vs uninf testing odds ratio

print(B)
print(Phi)

Y_0 <- 1e-4      ## initial prevalence
N <- 1e6         ## pop size
NY_0 <- N*Y_0    ## initial number infected

mp_default(sir_season)$matrix

# r <- log(2)/3    ## growth rate (doubling time = 3)
beta_high <- 0.25  ## Highest seasonal beta value
beta_low <- 0.05   ## Lowest seasonal beta value
gamma <- 0.1       ## Recovery rate
eta <- 0.02        ## Immunity waning rate
period <- 365/2    ## seasonal period
tmin <- 365*5      ## start when arrive at the stable phase
tmax <- 365*7      ## 2 years, 4 period observation 
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
# true_param <- tibble::lst(    log_B=log(B)
#                             , log_Phi=log(Phi)
#                             , logY_0=log(Y_0)
#                             , logit_T_B=logit_trans(T_B)
#                             , logit_T_Y=logit_trans(T_Y)
# )

tp_list <-tibble::lst(  beta_low
                      , beta_high
                      , period
                      , gamma
                      , eta
                      , N
                      , T_B
                      , T_Y
                      , I = NY_0
                      , R = 0
)

(sir_season
  |> mp_tmb_update(
    default = tp_list
  )
  |> mp_tmb_insert_backtrans(variables = c("beta_low","beta_high","gamma","eta","period","I"), mp_log)
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
)->sirs_seasonal

sirs_seasonal |> mp_expand()

(sirs_seasonal
  |> mp_simulator(
    time_steps = tmax
    , outputs = c("pY","T_prop","pos","OT","OP","I","S")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  
) -> dat_all

### only observe data from year 5 to year 7
dat<- (dat_all|> filter(time>=tmin)
)

###
(dat 
  |> filter(matrix=="S"|matrix=="I"|matrix=="OP"|matrix=="OT")
) -> dat_count

(dat 
  |> filter(matrix=="pY"|matrix=="pos"|matrix=="T_prop")
) -> dat_ratio


print(ggplot(dat_count)
      + aes(time, value, color=matrix)
      + geom_line()
      + scale_y_log10()
)

print(ggplot(dat_ratio)
      + aes(time, value, color=matrix)
      + geom_line()
      + scale_y_log10()
)

############################################

##### Start values of the fitting
### Assume we "know" S
## True $S$ value at t_min
S <- dat_all[which(dat_all$time==tmin-1 & dat_all$matrix=="S"),]$value

## hat_S as the starting point
hat_S <- S*(1)

print(S)
print(hat_S)

### approximation of hat{I} inferred from first data point:
OT <- dat[which(dat$time==tmin & dat$matrix=="OT"),]$value
OP <- dat[which(dat$time==tmin & dat$matrix=="OP"),]$value

print(T_Y)
### we don't know the true T_Y but have some estimation
hat_T_Y <- T_Y

### Inferred value from assumptions
hat_T <- OT/N
hat_p <- OP/OT
hat_Y <- (hat_p*hat_T)/hat_T_Y

hat_I <- hat_Y*N

### In fitting, should check if S+I >=1
if(hat_S+hat_I>N){
  print("initial S+I value larger than N")
} else {"check"}

### starting values for fitting
# print(tp_list)
sp_list <-tibble::lst(  beta_low=beta_low
                      , beta_high=beta_high
                      , period=period
                      , gamma=gamma
                      , eta=eta
                      , N=N
                      , T_B=T_B
                      , T_Y=hat_T_Y
                      , S=hat_S
                      , I=hat_I
                      )

### Change simulation model
(sir_season
  |> mp_tmb_update(
    default = sp_list
  )
  |> mp_tmb_insert_backtrans(variables = c("beta_low","beta_high","gamma","eta","period","S","I"), mp_log)
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
    , at = 10
    , expressions = list(
      R ~ N - S - I
    )
  )
  ## |> mp_tmb_delete(phase = "before", at = Inf, default = c("beta","gamma","I","T_B","T_Y"))
) -> sirs_seasonal_sim

print(sirs_seasonal_sim)

fit_pars <- c(  "log_beta_low"
              , "log_beta_high"
              , "log_gamma"
              , "log_eta"
              , "log_period"
              , "log_S"
              , "log_I"
              , "logit_T_B"
              , "logit_T_Y"
              )

calibrator <- mp_tmb_calibrator(
  sirs_seasonal_sim
  , data = dat
  , traj = c("OT", "OP", "T_prop", "pos")
  , par = fit_pars,
  , default = list(N=N)
)
## modify likelihood function (eventually we'll have mp_binom() so we can do this
## when defining the calibrator)
calibrator$simulator$replace$obj_fn(~ - sum(dbinom(obs_OT, N, sim_T_prop)) - sum(dbinom(obs_OP, obs_OT, sim_pos)))

calibrator|>print()

#mp_optimize(calibrator,optimizer = "optim", method = "BFGS")
#mp_optimize(calibrator)


fit<-mp_optimize(calibrator)
names(fit$par)<-fit_pars
print(fit)
exp(fit$par)
eta

fit_bk<-c( exp(fit$par[1])
          ,exp(fit$par[2])
          ,exp(fit$par[5])
          ,exp(fit$par[3])
          ,exp(fit$par[4])
          ,N=N
          ,logit_backtrans(fit$par[8])
          ,logit_backtrans(fit$par[9])
          ,exp(fit$par[6])
          ,exp(fit$par[7])
          )
names(fit_bk)<-names(sp_list)
# sp_list
# fit$par
fit_bklist<-as.list(append(fit_bk,fit$par))


fit_bklist
print(fit_bk)

(sirs_seasonal_sim
  |> mp_tmb_update(phase = "during", default = fit_bklist)
  |> mp_simulator(
    time_steps = tmax-tmin+1
    , outputs = c(  "OP"
                  , "OT"
                  , "pY"
                  , "S"
    )
  ) 
)->sirs_optim

sim_vals <-(sirs_optim|> mp_trajectory()
            |> dplyr::select(-c(row, col))
            # |> filter(time>=tmin)
)
# mp_default(sirs_optim)


sirs_optim
fit_bklist
sim_vals[1,]

# (dat_sim
#   |> pivot_wider(names_from = matrix,values_from = value)
# ) |> print(n=pts)
dat_tp<-filter(dat,matrix=="OP"|matrix=="OT"|matrix=="pY"|matrix=="S")

# dat_sim$time<-dat_sim$time+tmin-1
sim_vals$time<-sim_vals$time+tmin-1

# dat_sim <- cbind(dat_sim,model=rep("sim_init",length(dat_sim[,1])))
dat_tp <- cbind(dat_tp,model=rep("real",length(dat_tp[,1])))
sim_vals <- cbind(sim_vals,model=rep("optim",length(sim_vals[,1])))

dat_compare<-rbind(dat_tp,sim_vals)
fit_curve <- (ggplot(dat_compare)
              + aes(x=time, y=value, color=matrix, linetype = model)
              + geom_line()
              #+ geom_point(data=dat_tp,aes(shape="real"))
              + scale_y_log10()
              #+ scale_colour_manual(values = c("blue", "red", "black"))
              + scale_linetype_manual(values = c(1,2,3))
)
print(fit_curve)

dat_tp[1:5,]
sim_vals[1:5,]

fit_tp_result <- (mp_tmb_coef(calibrator,conf.int = TRUE) 
                  |> select(-c("term", "type","row","col","std.error"))
                  |> cbind(true_value=c(  beta_low = beta_low
                                        , beta_high = beta_high
                                        , gamma = gamma
                                        , eta = eta
                                        , period = period
                                        , S=S
                                        , I=dat_all[which(dat_all$time==tmin-1 & dat_all$matrix=="I"),]$value
                                        , T_B=T_B
                                        , T_Y=T_Y)
                  )
)
print(fit_tp_result)
fit_compare <- ggplot(fit_tp_result, aes(y = mat)) +
  geom_pointrange(aes(x = estimate, xmin = conf.low, xmax = conf.high)) +
  geom_point(aes(x=true_value), colour = "red") +
  facet_wrap(~mat, ncol = 1, scale  = "free")
print(fit_compare)

