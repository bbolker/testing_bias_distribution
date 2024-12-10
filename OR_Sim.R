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

## Simulate the data
dat <- tibble(t=t
	, pY = pmin(pY_0*exp(r*t), 1)
	, NY = rbinom(pts, N, pY)
	, posTests = rbinom(pts, NY, T_Y)
	, negTests = rbinom(pts, N-NY, T_B)
)

long_dat <- (dat
	|> select(-pY)
	|> pivot_longer(-t)
)

# print(ggplot(long_dat)
# 	+ aes(t, value, color=name)
# 	+ geom_line()
# 	+ scale_y_log10()
# )

### function to calculate negative log-likelihood:
LL <- function(B,Phi,pY_0,r,dat,N,tmax){
  T_B <- B/(1+B)
  T_Y <- B*Phi/(1+B*Phi)
  t <- c(0:tmax)
  ### simulated time series
  sim <- tibble(t=t
                , pY = pmin(pY_0*exp(r*t), 1)
                , NY = rbinom(pts, N, pY)
  )
  
  out <- (-sum(dbinom(dat$posTests, sim$NY, T_Y,log = TRUE))
          -sum(dbinom(dat$negTests,N-sim$NY,T_B,log = TRUE)))
  return(out)
}

real_ML<-LL(B,Phi,pY_0,r,dat,N,tmax)
real_ML

