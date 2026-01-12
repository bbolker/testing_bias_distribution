library(shellpipes)

seed <- 0376

N <- 1000    #population
I0 <- 1      #initial infection
beta <- 0.5  #transmission rate
D <- 6       #recovery duration
deltat <- 1  # discrete sim step size
# steps <- 100 # discrete sim total step

tmin <- 30
tmax <- 50  

## Test ratio ?? Relative hazard?
h <- 2

## awareness parameters
w0 <- 0.2
wI <- 5
alpha <- 3

## Backward compatibility; deprecated
rho <- 0.3
sig <- 0.05

## dispersion for starting points
cn <- 0.2

saveEnvironment()
