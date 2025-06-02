library(macpan2)

## flow diagram specification
flows = list(
    time ~ time+1
  , beta ~ (1/2)*(cos((time-1)*2*pi/(period))+1)*(beta_0) 
  , mp_per_capita_flow("S", "I", "beta * I / N", "infection")
  , mp_per_capita_flow("I", "R", "gamma", "recovery")
  , mp_per_capita_flow("R", "S", "phi", "waning_immunity")
)

## default values for quantities required to run simulations
default = list(
    beta_0 = 0.2    ## peak transmission rate for each season
  , period = 365/2  ## period of each season 
  , gamma = 0.05    ## recovery rate
  , pi = 3.141593   ## pi value for seasonal transmission
  , phi = 0.01      ## immunity wanning rate
  , N = 10000       ## total population size (constant in this model)
  , I = 1           ## initial number of infectious individuals
  , R = 0           ## initial number of recovered individuals
  , time = 0        ## time variable for the time varying beta
)

## compute the initial number of susceptible individuals
initialize_state = list(S ~ N - I - R)

## model specification
spec = mp_tmb_model_spec(
    before = initialize_state
  , during = flows
  , default = default
)
