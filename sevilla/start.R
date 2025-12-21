library(shellpipes)

loadEnvironments()

bstart <- rlnorm(1, log(beta), cn)
Dstart <- rlnorm(1, log(D), cn)
Istart <- rlnorm(1, log(I0), cn)
hstart <- rlnorm(1, log(h), cn)

saveEnvironment()
