library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)

# Set Seeds
set.seed(13521)

## Initial true values:
T_B <- 0.04
T_Y <- 0.5
B <- T_B/(1-T_B)
Phi <- (T_Y/(1-T_Y))/B

print(B)
print(Phi)

Y_0 <- 1e-4
N <- 1e6
NY_0 <- N*Y_0

r <- log(2)/3
tmax <- 19
t <- c(0:tmax)
pts <- length(t)

## Simulate the data
dat <- tibble(t=t
	, Y = pmin(Y_0*exp(r*t), 1)
	## , NY = rbinom(pts, N, Y)
	, NY = round(N*Y)
	, posTests = rbinom(pts, NY, T_Y)
	, negTests = rbinom(pts, N-NY, T_B)
)

long_dat <- (dat
	|> select(-Y)
	|> pivot_longer(-t)
)

# print(ggplot(long_dat)
# 	+ aes(t, value, color=name)
# 	+ geom_line()
# 	+ scale_y_log10()
# )

### function to calculate negative log-likelihood:
nll <- function(B,Phi,logY_0,r,dat,N,tmax){
	Y_0 <- exp(logY_0)
  T_B <- B/(1+B)
  T_Y <- B*Phi/(1+B*Phi)
  t <- c(0:tmax)
  pts <- length(t)
  Ynull = Y_0*exp(r*t) ## Going to truncate here for now; should we change sim instead?
  ### simulated time series
  sim <- tibble(t=t
                # , Y = 1-exp(-Ynull)
                , Y = pmin(Ynull, 1)
                , NY = round(N*Y)
                # , NY = rbinom(pts, N, Y)
  )
  if (max(sim$NY<dat$posTests)) cat("Positive tests exceed infected population", "\n")
  if (max((N-sim$NY)<dat$negTests)) cat("Negative tests exceed uninfected population", "\n")
  out <- (-sum(dbinom(dat$posTests, sim$NY, T_Y,log = TRUE))
          -sum(dbinom(dat$negTests,N-sim$NY,T_B,log = TRUE)))
  return(out)
}

real_ML<-nll(B,Phi,log(Y_0),r,dat,N,tmax)
real_ML

print(nll(B,Phi,log(0.5),0.0,dat,N,tmax))

mle_out_debug <- try(mle2(nll
                          ,start = list(B=B
                                        ,Phi=Phi
                                        ,logY_0=log(Y_0)
                                        ,r=r)
                          ,data = list(dat=dat
                                       ,N=N
                                       ,tmax=tmax)
                          ,control = list(maxit=10000)
                          ,method = "Nelder-Mead"
))
warnings()

print(mle_out_debug)
