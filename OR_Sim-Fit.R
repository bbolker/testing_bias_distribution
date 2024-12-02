library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)

set.seed(1351)

## Initial true values:
T_B <- 0.04
T_Y <- 0.5
B <- T_B/(1-T_B)
Phi <- (T_Y/(1-T_Y))/B

print(B)
print(Phi)

NY_0 <- 100
N <- 1e6
r <- log(2)/3

true_pars <- c(NY_0=NY_0, r=r, B=B, Phi=Phi, N=1e6, T_B=T_B, T_Y=T_Y)
true_pars
## simulation
maxtime <- 29
t<-c(0:maxtime)

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
LL <- function(B,Phi,Y0,r,dd,ture_pars,tmax){
  T_B <- B/(1+B)
  T_Y <- B*Phi/(1+B*Phi)
  t <- c(0:tmax)
  N <- true_pars["N"]
  
  NY_t <- N*Y0*exp(r*t)
  # prevalence
  Y_t <- NY_t/N
  
  # Data frame to storage the data
  df<- data.frame(t,NY_t,Y_t)
  
  # expected test proportion T
  df$T <- (1-df$Y_t)*T_B+df$Y_t*T_Y
  
  ### number of test as parameter:
  df$TNum <- round(df$T*N)
  
  # expected test positivity P
  df$P <- df$Y_t*T_Y/df$T
  
  -sum(stats::dpois())
}

