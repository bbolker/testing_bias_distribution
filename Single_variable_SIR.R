library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(pracma)
library(gsl)
library(deSolve)
source("testing_funs.R")

#### An SIR trajectory simulator
MASIR_Proc <- function(beta, gamma, init_S=1e-3, ODEmaxTime=50, ODEstep=1e-2,TrackDyn=T){
  if (TrackDyn==T){
    Sys <- function(t, y, parms){
      with(as.list(c(parms,y)),{
        dS <- (-b)*I*S
        dI <- b*I*S-g*I
        dR <- g*I
        return(list(c(dS,dI,dR)))
      })
    }
    parms <- c(b=beta,g=gamma)
    times <- seq(0,ODEmaxTime,by=ODEstep)
    y <- c(S=init_S,I=1-init_S,R=0)
    
    Sys_out <- ode(y,times,Sys,parms)
  }
  
  g <- gamma
  b <- beta
  S_0 <- init_S
  R_0 <- 0
  
  R0 <- b/g
  RInf <- Sys_out[length(Sys_out[,4]),4]
  
  if (TrackDyn==T){
    return(list(R0=R0,RInfinity=RInf, Dynamic=Sys_out))
  } else {
    return(list(R0=R0,RInfinity=RInf))
  }
}

maxtime<-150
ODEstep<-1e-1

SIR_Sys<- MASIR_Proc(0.40, 0.20, init_S = 0.999, ODEmaxTime=maxtime, ODEstep=ODEstep,TrackDyn = T)
SIR_Sys$R0
SIR_Sys$RInfinity
SIR_Curve <- SIR_Sys$Dynamic

timeseq <- seq(1,(maxtime+1)*(1/ODEstep),(1/ODEstep))
Curve_tseq <- SIR_Curve[timeseq,]
I_tseq <- Curve_tseq[,3]
# SIR_inc <- rep(0,length(timeseq)-1)
# #plot(Sys_Opt[,1],Sys_Opt[,3])
# 
# for (i in c(1:(length(timeseq)-1))){
#   SIR_inc[i] <- Curve_tseq[i+1,4]-Curve_tseq[i,4]
# }

phi_const <- 0.2
phi_seq <- rep(phi_const,length(timeseq))
test_seq<-runif(length(timeseq),min=0.001,max = 0.01)

prop_seq <- prop_pos_test_new(I_tseq,test_seq,phi_seq,method="log")

prop_seq
plot(c(1:(maxtime+1)),prop_seq,type="l")


### If knowing phi, we area able to figure out inc

######### Single variable example
inc_fn <- function(inc,phi,test_obs,prop_obs){
  prop_pred <- prop_pos_test_new(inc,test_obs,phi,method="log")
  out <- prop_pred-prop_obs
  return(out)
}

#E.g. t=25
tp <- 25
inc_root <- uniroot(inc_fn,c(0.00001,0.5),phi=phi_seq[tp],test_obs=test_seq[tp],prop_obs=prop_seq[tp])
inc_root$root
I_tseq[tp]

# Reasonable time points
inc_root_seq <- c(1:100)
for(i in c(1:100)){
  out <- uniroot(inc_fn,c(0.00001,0.5),phi=phi_seq[i],test_obs=test_seq[i],prop_obs=prop_seq[i])
  inc_root_seq[i] <- out$root
}

plot(c(1:100),I_tseq[1:100],type="l")
points(c(1:100),inc_root_seq,col="red")

### vise versa for phi
phi_fn <- function(phi,inc,test_obs,prop_obs){
  prop_pred <- prop_pos_test_new(inc,test_obs,phi,method="log")
  out <- prop_pred-prop_obs
  return(out)
}

#E.g. t=25
tp <- 25
phi_root <- uniroot(phi_fn,c(0.01,0.5),inc=I_tseq[tp],test_obs=test_seq[tp],prop_obs=prop_seq[tp])
phi_root$root

# Reasonable time points: i not too small
phi_root_seq <- c(1:100)
for(i in c(1:100)){
  out <- uniroot(phi_fn,c(0.01,0.5),inc=I_tseq[i],test_obs=test_seq[i],prop_obs=prop_seq[i])
  phi_root_seq[i] <- out$root
}
plot(c(1:100),phi_seq[1:100],type="l")
points(c(1:100),phi_root_seq,col="red")

