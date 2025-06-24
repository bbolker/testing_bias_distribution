library(macpan2)

## flow diagram specification
flows = list(
    beta ~ beta_low+(1/2)*(cos((time_step(0)+phase)*2*pi/(period))+1)*(beta_high-beta_low) 
  , mp_per_capita_flow("S", "I", "beta * I / N", "infection")
  , mp_per_capita_flow("I", "R", "gamma", "recovery")
  , mp_per_capita_flow("R", "S", "eta", "waning_immunity")
)

## default values for quantities required to run simulations
default = list(
    beta_high = 0.2 ## Highest transmission rate for each season
  , beta_low = 0.02 ## Lowest transmission rate for each season
  , period = 365/2  ## period of each season 
  , gamma = 0.05    ## recovery rate
  , pi = pi         ## pi value for seasonal transmission
  , eta = 0.01      ## immunity waning rate
  , N = 10000       ## total population size (constant in this model)
  , I = 1           ## initial number of infectious individuals
  , R = 0           ## initial number of recovered individuals
  , phase = 0       ## the phase t of time varying beta(t) at the current time
                    ## can only be positive integer
)

## compute the initial number of susceptible individuals
initialize_state = list(S ~ N - I - R)

## model specification
spec = mp_tmb_model_spec(
    before = initialize_state
  , during = flows
  , default = default
)
