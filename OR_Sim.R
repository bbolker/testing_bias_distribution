library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)

# Set Seeds
set.seed(13521)

## Initial true values:
T_B <- 0.04
T_Y <- 0.5
B <- T_B/(1-T_B)
Phi <- (T_Y/(1-T_Y))/B

print(B)
print(Phi)

pY_0 <- 1e-4
N <- 1e6
r <- log(2)/3
tmax <- 29
t <- c(0:tmax)
pts <- length(t)

dat <- tibble(t=t
	, pY = pmin(pY_0*exp(r*t), 1)
	, NY = rbinom(pts, N, pY)
	, posTests = rbinom(pts, NY, T_Y)
	, negTests = rbinom(pts, N-NY, T_B)
)

long <- (dat
	|> select(-pY)
	|> pivot_longer(-t)
)

print(ggplot(long)
	+ aes(t, value, color=name)
	+ geom_line()
	+ scale_y_log10()
)
