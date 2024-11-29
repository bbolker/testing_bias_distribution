library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)

set.seed(1351)

## Initial true values:
T_B <- 0.2
T_Y <- 0.5
B <- T_B/(1-T_B)
Phi <- (T_Y/(1-T_Y))/B

print(B)
print(Phi)

Y_0 <- 10
N <- 1e6
r <- log(2)/3

## simulation
t<-c(0:29)
NY_t <- Y_0*exp(r*t)


