library(shellpipes)

seed <- 0376

N <- 1000
I0 <- 1
beta <- 0.2
D <- 10
deltat <- 1
steps <- 100

## Test ratio
h <- 2

## awareness parameters
w0 <- 0.
wI <- 5
alpha <- 3

## Backward compatibility; deprecated
rho <- 0.3
sig <- 0.05

## dispersion for starting points
cn <- 0.2

saveEnvironment()
