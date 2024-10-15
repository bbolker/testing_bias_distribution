library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(DPQ)
library(deSolve)
source("testing_funs.R")

### Numeric Check for explicit soln of p(t)
a<-1
b<-seq(from=40,to=50,by=0.1)
t<-0.15
# t <- 0.15
##t <- 0.25

# quick testing of simplified expression
int_p<-rep(0,length(b))
cdf_p<-rep(0,length(b))
simp_p <- rep(0,length(b))
log_simp_p <- rep(0,length(b))
for (i in c(1:length(b))) {
  lwr <- Qbeta2(t, a, b[i],lower.tail=FALSE)
  int_p[i] <- (integrate(\(x) dbeta(x, a, b[i])*x, lwr, 1)$value)/t
  cdf_p[i] <- (pbeta(lwr,a+1,b[i],lower.tail=FALSE)*(a/(a+b[i])))/t
  simp_p[i] <- (a/(a+b[i])*(t+(lwr^a*(1-lwr)^b[i])/(beta(a,b[i])*a)))/t
  log_simp_p[i] <- (a/(a+b[i])*(t+(exp(a*log(lwr)+b[i]*log(1-lwr)))/(beta(a,b[i])*a)))/t
}

# plot(b,cdf_p)
# points(b,int_p,col="blue",pch=6)
# points(b,simp_p,col="red",pch=3)
# points(b,log_simp_p,col="green",pch=3)


### log.space difference between different methods
## Idea: heat map on inc \in (0,1) & phi \in (0,1)
## ? how to deal with negative difference
## Current solution, just consider absolute value of difference?
## Testing ideas
out <- c(1:length(b))
for (i in c(1:length(b))) {
  a <- sign (log(simp_p[i])-log(log_simp_p[i]))
  if (a>=0) {
    out[i] <- logspace.sub(log(simp_p[i]),log(log_simp_p[i]))
  } else {
    out[i] <- logspace.sub(log(log_simp_p[i]),log(simp_p)[i])
  }
}

## define function abs_log_sub:
abs_log_sub <- function(x,y,zeroHL=-Inf){
  sgn <- sign(log(x)-log(y))
  ## Exclude failing case: possible reason, Qbeta gives 1 and Pbeta generate infinity
  ##### FIXIT??? 
  if (is.infinite(sgn)){
    out <- NaN
    } else if(is.nan(sgn)){
    out <- NaN
    } else {
    if (sgn>0) {
      out <- logspace.sub(log(x),log(y))
    } else if (sgn==0) {
      out <- logspace.sub(log(x),log(y))
      if (out==-Inf){
        out <- zeroHL
      }
      # Highlight very close to (0) results
    } else {
      out <- logspace.sub(log(y),log(x))
    }
  }
  return(out)
}
abs_log_sub <- Vectorize(abs_log_sub,c("x","y"))

abs_log_sub(simp_p,log_simp_p)
## Heat map generating
## integration "int" method seems failed for some extreme values: generate infinite values?
dd <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                   phi=seq(from=0.01,to=0.99,by=0.01),
                   inc=seq(from=0.01,to=0.99,by=0.01))
       %>% as_tibble()
       #%>% mutate(pos_prop_int=prop_pos_test_new(inc,test_prop,phi,method="int"))
       %>% mutate(pos_prop_cdf=prop_pos_test_new(inc,test_prop,phi,method="cdf"))
       %>% mutate(pos_prop_simp=prop_pos_test_new(inc,test_prop,phi,method="simp"))
       %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
       #%>% mutate(diff_cdf_log=sign(pos_prop_cdf-pos_prop_log))
       %>% mutate(diff_cdf_log=abs_log_sub(pos_prop_cdf,pos_prop_log,zeroHL=-45))
       %>% mutate(diff_simp_log=abs_log_sub(pos_prop_simp,pos_prop_log,zeroHL=-40))
       )

print(ggplot(dd,aes(phi,inc,fill=diff_cdf_log,group=test_prop))
      + geom_tile()
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_fill_viridis_c(expand=c(0,0))
      + ggtitle("log-space diff between cdf and log_simp")
)

print(ggplot(dd,aes(phi,inc,fill=diff_simp_log,group=test_prop))
      + geom_tile()
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_fill_viridis_c(expand=c(0,0))
      + ggtitle("log-space diff between simp and log_simp")
)

#### need to investigate more: what happens to a phi=0.99 case  
prop_pos_test_new(0.5,0.001,0.97,method = "cdf",debug = T)
prop_pos_test_new(0.5,0.001,0.97,method = "log",debug = T)
prop_pos_test_new(0.5,0.001,0.97,method = "simp",debug = T)

val_log <- function(i,phi,t) {
  phi_0 <- phi
  phi <- -log(1-phi)
  a <- i/phi; b <- (1-i)/phi
  qb<-qbeta(t,a,b, lower.tail = FALSE)
  log_qb <- log(qb)
  simp<-a/(a+b)*(t+(exp(a*log(qb)+b*log(1-qb)))/(beta(a,b)*a))
  log_simp <- log(simp)
  cdf <- pbeta(qb,a+1,b,lower.tail=FALSE)*(a/(a+b))
  log_cdf <- log(cdf)
  return(c(qb,log_qb,log_simp,log_cdf))
}

val_log(0.3,0.99,0.001)
val_log(0.3,0.98,0.001)
val_log(0.3,0.97,0.001)

####### Numerical observation
## even lwr output is "1", log(qbeta) could still have difference
## using both simp and log simp could somehow catch the difference
## However, for phi=0.99 very extreme case, even log(qbeta) gets solid zero

## very extreme case have positivity=i，since the log(1-qb)=-Inf
## cdf method will give "solid" zero for that case (log=-Inf)