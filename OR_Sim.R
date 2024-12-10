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

pY_0 <- 1e-4
N <- 1e6
NY_0 <- N*pY_0

r <- log(2)/3
tmax <- 19
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
  pts <- length(t)
  ### simulated time series
  sim <- tibble(t=t
                , pY = pmin(pY_0*exp(r*t), 1)
                ### round here???
                , NY = round(N*pY)
                # , NY = rbinom(pts, N, pY)
  )
  if (max(sim$NY<dat$posTests)) cat("Underestimated infected population, positive tests exceed infected population", "\n")
  if (max((N-sim$NY)<dat$negTests)) cat("Overestimated infected population, negative tests exceed uninfected population", "\n")
  out <- (-sum(dbinom(dat$posTests, sim$NY, T_Y,log = TRUE))
          -sum(dbinom(dat$negTests,N-sim$NY,T_B,log = TRUE)))
  return(out)
}

real_ML<-LL(B,Phi,pY_0,r,dat,N,tmax)
real_ML

LL(0.05,Phi,pY_0,0.2,dat,N,tmax)

mle_out_debug <- try(mle2(LL
                          ,start = list(B=0.06
                                        ,Phi=Phi
                                        ,pY_0=pY_0
                                        ,r=0.2)
                          ,data = list(dat=dat
                                       ,N=N
                                       ,tmax=tmax)
                          ,control = list(maxit=10000)
                          ,method = "Nelder-Mead"
))
warnings()
r
mle_out_debug
# mle_out_NM <- update(mle_out, method = "Nelder-Mead")

######### The optimization always leads to negative pY_0

mle_out2 <- try(mle2(LL
                          ,start = list(B=B
                                        ,Phi=50
                                        ,pY_0=9e-5
                                        ,r=r)
                          ,data = list(dat=dat
                                       ,N=N
                                       ,tmax=tmax)
                          ,control = list(maxit=10000)
                          ,method = "Nelder-Mead"
))

mle_out2

test_params <- list(B=B
                    ,Phi=22
                    ,pY_0=8e-5
                    ,r=r)

mledat <- list(dat=dat
               ,N=N
               ,tmax=tmax)

do.call(LL, c(test_params, mledat))


### Heatmap?

