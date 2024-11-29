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

Y_0 <- 100
N <- 1e6

