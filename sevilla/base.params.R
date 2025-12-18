library(shellpipes)

seed <- 0376

N <- 1000
I0 <- 1
beta <- 0.2
D <- 10
deltat <- 1
steps <- 100

## Test ratio
h <- 5

## awareness parameters
w0 <- 0.1
wI <- 1
alpha <- 1

## Backward compatibility; deprecated
rho <- 0.3
sig <- 0.05

saveEnvironment()
