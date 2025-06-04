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
                   , default=list(gamma=0.02) 
                   )
  |> mp_simulator(
      time_steps = 2000
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
                   , default=list(  gamma=0.02
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

