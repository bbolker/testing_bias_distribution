library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(pracma)
library(gsl)
library(deSolve)
source("testing_funs.R")
source("Single_variable_SIR.R")

maxtime<-150
ODEstep<-1e-1

SIR_Sys<- MASIR_Proc(0.40, 0.20, init_S = 0.999, ODEmaxTime=maxtime, ODEstep=ODEstep,TrackDyn = T)
SIR_Sys$R0
SIR_Sys$RInfinity
SIR_Curve <- SIR_Sys$Dynamic

timeseq <- seq(1,(maxtime+1)*(1/ODEstep),(1/ODEstep))
Curve_tseq <- SIR_Curve[timeseq,]
I_tseq <- Curve_tseq[,3]

phi_const <- 0.2
phi_seq <- rep(phi_const,length(timeseq))
set.seed(1746)
test_seq<-runif(length(timeseq),min=0.001,max = 0.01)

prop_seq <- prop_pos_test_new(I_tseq,test_seq,phi_seq,method="log")

prop_seq
plot(c(1:(maxtime+1)),prop_seq,type="l")

### LSS optimize
target_fn <- function(x,test_obs,prop_obs){
  inc_val <- x[1]
  phi_val <- x[2]
  ## LSS
  prop_pred1 <- prop_pos_test_new(inc_val,test_obs[1],phi_val,method="log")
  prop_pred2 <- prop_pos_test_new(inc_val,test_obs[2],phi_val,method="log")
  LSS <- (prop_pred1-prop_obs[1])^2+(prop_pred2-prop_obs[2])^2
  return(LSS)
}

test_case <- nlm(target_fn,c(0.01,0.2),test_obs=c(test_seq[55],test_seq[56]),prop_obs=c(prop_seq[55],prop_seq[56]), stepmax=0.008,iterlim = 1000,gradtol = 1e-12)

test_case$estimate

I_tseq[55]
I_tseq[56]

####### Fail to identify both phi and inc

# a<-prop_pos_test_new(test_case$estimate[1], test_seq[1], test_case$estimate[2], method = "log")
# b<-prop_pos_test_new(test_case$estimate[1], test_seq[2], test_case$estimate[2], method = "log")
# 
# c<-prop_pos_test_new(I_tseq[1], test_seq[1], phi_const, method = "log")
# d<-prop_pos_test_new(I_tseq[2], test_seq[2], phi_const, method = "log")
# 
# (a-c)^2+(b-d)^2

#### Brute-force heat map
target_fn2 <- function(inc,phi,test_obs,prop_obs){
  inc_val <- inc
  phi_val <- phi
  ## LSS
  prop_pred1 <- prop_pos_test_new(inc_val,test_obs[1],phi_val,method="log")
  prop_pred2 <- prop_pos_test_new(inc_val,test_obs[2],phi_val,method="log")
  LSS <- (prop_pred1-prop_obs[1])^2+(prop_pred2-prop_obs[2])^2
  return(LSS)
}

dd <- (expand.grid(phi=seq(from=0.001,to=0.5,by=0.001),
                   inc=seq(from=0.0001,to=0.05,by=0.0001))
       %>% as_tibble()
       %>% mutate(LSS=target_fn2(inc,phi,test_obs=c(test_seq[55],test_seq[56]),prop_obs=c(prop_seq[55],prop_seq[56])))
       )

print(ggplot(dd,aes(phi,inc,fill=log(LSS)))
      + geom_tile()
      #+ scale_y_log10()
      + scale_fill_viridis_c(expand=c(0,0))
      + ggtitle("LSS of phi and inc")
)

dd[which.min(dd$LSS),]
test_case$estimate

I_tseq[55]
I_tseq[56]

## If knowing phi=0.2
Opt_knowphi <-which.min(dd$LSS[which(dd$phi==phi_const)])
dd[which(dd$phi==phi_const)[Opt_knowphi],]

I_tseq[55]
I_tseq[56]

### If knowing phi, we area able to figure out inc

######## Time series:
## Idea: take 2-3 days as a whole, test_prop1=day1 test_prop2=day1+day2, test_prop3=day1+day2+day3
## data: tp1,tp2,tp3, pp1,pp2,pp3
## Input: test_prop1= tp1, test_prop2= tp1+tp2, test_prop3=tp1+tp2+tp3
## pos_prop1=pp1, pos_prop2= (tp1*pp1+tp2*pp2)/(tp1+tp2), test_prop3=(tp1*pp1+tp2*pp2+tp3*pp3)/(tp1+tp2+tp3)
## trying to figure out i (and phi) 

## Need to think more carefully of the assumption
#### What is the impact of time theories/aggragation

target_fn3 <- function(inc,phi,test_obs,prop_obs){
  inc_val <- inc
  phi_val <- phi
  
  test_prop1 <- test_obs[1]
  test_prop2 <- test_obs[1]+test_obs[2]
  test_prop3 <- test_obs[1]+test_obs[2]+test_obs[3]
  
  pos_prop1 <- prop_obs[1]
  pos_prop2 <- (prop_obs[1]*test_obs[1]+prop_obs[2]*test_obs[2])/test_prop2
  pos_prop3 <- (prop_obs[1]*test_obs[1]+prop_obs[2]*test_obs[2]+prop_obs[3]*test_obs[3])/test_prop3
  
  ## LSS
  prop_pred1 <- prop_pos_test_new(inc_val,test_prop1,phi_val,method="log")
  prop_pred2 <- prop_pos_test_new(inc_val,test_prop2,phi_val,method="log")
  prop_pred3 <- prop_pos_test_new(inc_val,test_prop3,phi_val,method="log")
  
  LSS <- (prop_pred1-pos_prop1)^2+(prop_pred2-pos_prop2)^2+(prop_pred3-pos_prop3)^2
  return(LSS)
}


time <- 15
test_obs=c(test_seq[time],test_seq[time+1],test_seq[time+2])
prob_obs=c(prop_seq[time],prop_seq[time+1],prop_seq[time+2])


dd2 <- (expand.grid(phi=seq(from=0.001,to=0.5,by=0.001),
                   inc=seq(from=0.001,to=0.1,by=0.001))
       %>% as_tibble()
       %>% mutate(LSS=target_fn3(inc,phi,test_obs=test_obs,prop_obs=prob_obs))
)

print(ggplot(dd2,aes(phi,inc,fill=log(LSS)))
      + geom_tile()
      #+ scale_y_log10()
      + scale_fill_viridis_c(expand=c(0,0))
      + ggtitle("LSS of phi and inc")
)

dd2[which.min(dd2$LSS),]

print(dd2[which(dd2$phi==0.2),],n=100)

I_tseq[time]
I_tseq[time+1]
I_tseq[time+2]


## not working.

