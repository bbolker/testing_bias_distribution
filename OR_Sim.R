
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
	, pY = pmin(Y_0*exp(r*t), 1)
	## , NY = rbinom(pts, N, pY)
	, NY = round(N*pY)
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
LL <- function(B,Phi,logY_0,r,dat,N,tmax, debug = TRUE){
	Y_0 <- exp(logY_0)
  T_B <- B/(1+B)
  T_Y <- B*Phi/(1+B*Phi)
  t <- c(0:tmax)
  pts <- length(t)
  ### simulated time series
  sim <- tibble(t=t
                , pY = pmin(Y_0*exp(r*t), 1)
                ### round here???
                , NY = round(N*pY)
                # , NY = rbinom(pts, N, pY)
  )
  if (max(sim$NY<dat$posTests)) cat("Underestimated infected population, positive tests exceed infected population", "\n")
  if (max((N-sim$NY)<dat$negTests)) cat("Overestimated infected population, negative tests exceed uninfected population", "\n")
  postest_nll <- -sum(dbinom(dat$posTests, sim$NY, T_Y,log = TRUE))
  negtest_nll <- -sum(dbinom(dat$negTests,N-sim$NY,T_B,log = TRUE))
  out <- postest_nll + negtest_nll
  if (debug) cat(B, Phi, logY_0, r, postest_nll, negtest_nll,
                 out, "\n")
  return(out)
}

real_ML<-LL(B,Phi,log(Y_0),r,dat,N,tmax)
real_ML

LL(0.05,Phi,log(Y_0),0.2,dat,N,tmax)

mle_out_debug <- try(
    
    mle2(LL
        ,start = list(B=B
                     ,Phi=Phi
                     ,logY_0=log(Y_0)
                     ,r=r)
        ,data = list(dat=dat
                    ,N=N
                    ,tmax=tmax)
        ,control = list(maxit=10000, parscale = c(B, Phi, log(Y_0), r))
       , method = "Nelder-Mead"
         )

)

## final value:
## 0.04169482 23.91323 -9.2344 0.232682 76.43475 136.2932 212.7279

warnings()
print(r)
mle_out_debug
# mle_out_NM <- update(mle_out, method = "Nelder-Mead")

