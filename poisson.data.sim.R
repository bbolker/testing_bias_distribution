# Sys.setenv(LANG = "en")

library(shellpipes)
rpcall("const.data.Rout float.data.sim.R const.params.rda")
rpcall("float.data.Rout float.data.sim.R float.params.rda")
suppressPackageStartupMessages(library(dplyr))
library(tidyr)

suppressPackageStartupMessages(library(macpan2))
options(macpan2_verbose = FALSE)

#library(ggplot2); theme_set(theme_bw())
loadEnvironments()

#source("mle2_tidy.R")
#source("spec_trans_par.R")
#source("float.params.R")

# Set Seeds
set.seed(seed)

gamma <- deltat/D
t <- c(tmin:tmax)
pts <- length(t) ## number of time points

logit_trans <- function(x){
  log(x)-log(1-x)
}
logit_backtrans <- function(x){
  (1/(1 + exp(-x)))
}

## BMB: use self-naming list from tibble pkg
true_param <- tibble::lst(  log_h=log(h)
                          , log_I0=log(I0)
                          , log_W0=log(w0)
                          , log_WI=log(wI)
                          , log_alpha=log(alpha)
                          )
print(true_param)

tp_list <-tibble::lst(  beta
                      , gamma
                      , N
                      , h
                      , w0
                      , wI
                      , alpha
                      , I = I0
                      , R = 0
)

### import SIR from macpan-model-lib
mc_sir <- mp_tmb_library("starter_models","sir", package = "macpan2")

(mc_sir
  |> mp_tmb_update(
    default = tp_list
    )
  |> mp_tmb_insert_backtrans(  variables = c( "beta"
                                             ,"gamma"
                                             ,"I","h"
                                             ,"w0"
                                             ,"wI"
                                             ,"alpha")
                             , mp_log)
  |> mp_tmb_insert(
      phase = "during"
    , at = Inf
    , expressions = list(
        pY ~ I/N                          ## Prevalence based on SIR
      , pSus ~ S/N                        ## Susceptible proportion  
      , b ~ w0+wI*pY*exp(-alpha*(1-pSus)) ## Concern floating baseline hazard
      , T_B ~ 1-exp(-b)
      , T_Y ~ 1-exp(-b-h)
      # , T_prop ~ (1-pY)*T_B+pY*T_Y        ## Expected test proportion
      # , T_pos ~ pY*T_Y/T_prop             ## Expected test positivity
      , ONeg ~ rpois(T_B*(N-I))
      , OPos ~ rpois(T_Y*I)
    )
  )
  ## |> mp_tmb_delete(phase = "before", at = Inf, default = c("beta","gamma","I","T_B","T_Y"))
)->sir

# sir |> mp_expand()

(sir
  |> mp_simulator(
      time_steps = tmax
    # , outputs = c("pY","T_prop","T_pos","OT","OP","I","S","T_B")
    , outputs = c("pY","ONeg","OPos","I","S","T_B")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  
) -> dat_all

dat<- (dat_all|> filter(time>=tmin)
)

(dat
  |> pivot_wider(names_from = matrix,values_from = value)
) -> dat_fit

(dat_all
  |> pivot_wider(names_from = matrix,values_from = value)
) -> dat_pall

dat_pall
saveEnvironment()