library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)

# Set Seeds
set.seed(13521)

## Initial true values:
T_B <- 0.04
T_Y <- 0.5
B <- 1-T_B
Phi <- (1-T_Y)/B

print(B)
print(Phi)

NY_0 <- 100
N <- 1e6
r <- log(2)/3

true_pars <- c(NY_0=NY_0, r=r, B=B, Phi=Phi, N=1e6, T_B=T_B, T_Y=T_Y)
true_pars
## simulation
tmax <- 29
t<-c(0:tmax)

# Real number of infected
NY_t <- NY_0*exp(r*t)

# Real prevalence
Y_t <- NY_t/N

# Data frame to storage the data
dd<- data.frame(t,NY_t,Y_t)

# Real expected test proportion T
dd$T <- (1-dd$Y_t)*true_pars["T_B"]+dd$Y_t*true_pars["T_Y"]

# Real expected test positivity P
dd$P <- dd$Y_t*true_pars["T_Y"]/dd$T

# dd$P-(dd$Y_t*true_pars["Phi"]*(1+true_pars["B"])/(1+(true_pars["B"]+dd$Y_t)*true_pars["Phi"]-dd$Y_t))

### Observed testing number N T*:
dd$OTNum <-  with(c(as.list(true_pars), dd),
               rpois(length(t), round(T*N))
               )

### Observed positive count N T P*:
dd$OPNum <- with(c(as.list(true_pars), dd)
              ,rbinom(length(t), size = OTNum, prob=P)
              )
### Observed testing proportion
dd$OT <- dd$OTNum/N

### Observed testing positivity
dd$OP <- dd$OPNum/dd$OTNum

# dd

### function to calculate negative log-likelihood:
LL <- function(B,Phi,NY_0,r,dd,N,tmax){
  T_B <- 1-B
  T_Y <- 1-B*Phi
  t <- c(0:tmax)
  
  NY_t <- NY_0*exp(r*t)
  # prevalence
  Y_t <- NY_t/N
  
  # Data frame to storage the data
  df<- data.frame(t,NY_t,Y_t)
  
  # expected test proportion T
  df$T <- (1-df$Y_t)*T_B+df$Y_t*T_Y
  
  ### number of test as parameter:
  df$TNum <- round(df$T*N)
  
  # expected test positivity P
  df$P <- (df$Y_t*T_Y)/df$T
  
  out <- -sum(stats::dpois(dd$OTNum, df$TNum,log = TRUE))-sum(stats::dbinom(dd$OPNum,df$TNum,df$P,log = TRUE))
  return(out)
}

real_ML<-LL(B,Phi,NY_0,r,dd,N,tmax)

# mle_out <- mle2(LL
#      ,start = list(B=true_pars["B"]
#                    ,Phi=true_pars["Phi"]
#                    ,NY_0=true_pars["NY_0"]
#                    ,r=true_pars["r"])
#      ,data = list(dd=dd
#                   ,N=true_pars["N"]
#                   ,tmax=tmax)
#      ,control = list(maxit=1000)
#      )


mle_out <- mle2(LL
     ,start = list(B=0.2
                   ,Phi=0.7
                   ,NY_0=150
                   ,r=0.3)
     ,data = list(dd=dd
                  ,N=true_pars["N"]
                  ,tmax=tmax)
     ,control = list(maxit=1000)
     )
# mle_out
true_pars[c("B","Phi","NY_0","r")]
coef(mle_out)

mle_out@details$value
real_ML
# summary(mle_out)
