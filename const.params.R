library(shellpipes)
loadEnvironments()
seed <- 0376

N <- 100000    #population
I0 <- 1      #initial infection
beta <- 0.5  #transmission rate
D <- 6       #recovery duration
deltat <- 1  # discrete sim step size

tmin <- 20
tmax <- 50  

## Relative hazard
h <- 0.5

## awareness parameters
w0 <- 0.05
wI <- 0 #const baseline
#w0 <- 0.05
#wI <- 1 # time-varying baseline
alpha <- 3

## Backward compatibility; deprecated
rho <- 0.3
sig <- 0.05

## dispersion for starting points
cn <- 0.2

saveEnvironment()
