library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(DPQ)
library(deSolve)
library(Rmpfr)
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
dd <- (expand.grid(test_prop=c(0.001,0.01,0.1,0.25),
                   phi=seq(from=0.01,to=0.99,by=0.005),
                   inc=seq(from=0.01,to=0.99,by=0.005))
       %>% as_tibble()
       #%>% mutate(pos_prop_int=prop_pos_test_new(inc,test_prop,phi,method="int"))
       %>% mutate(pos_prop_cdf=prop_pos_test_new(inc,test_prop,phi,method="cdf"))
       %>% mutate(pos_prop_simp=prop_pos_test_new(inc,test_prop,phi,method="simp"))
       %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
       %>% mutate(pos_prop_est=prop_pos_test_new(inc,test_prop,phi,method="est"))
       #%>% mutate(diff_cdf_log=sign(pos_prop_cdf-pos_prop_log))
       %>% mutate(diff_cdf_log=abs_log_sub(pos_prop_cdf,pos_prop_log,zeroHL=-45))
       %>% mutate(diff_simp_log=abs_log_sub(pos_prop_simp,pos_prop_log,zeroHL=-40))
       %>% mutate(diff_est_log=abs_log_sub(pos_prop_log,pos_prop_est))
       )

# fig_logdiff_cdf_log <- (
#   ggplot(dd,aes(phi,inc,fill=diff_cdf_log,group=test_prop))
#   + geom_raster()
#   + facet_wrap(~test_prop,scale="free",labeller = label_both)
#   # + scale_y_log10()
#   + scale_fill_viridis_c(expand=c(0,0))
#   + ggtitle("log-space diff between cdf and log_simp")
# )
#fig_logdiff_cdf_log
#ggsave("log-diff_cdf-log_simp.png",plot=fig_logdiff_cdf_log, path = "./pix", width=3200,height=1800,units="px")


fig_logdiff_simp_log <-(
  ggplot(dd,aes(phi,inc,fill=diff_simp_log,group=test_prop))
  + geom_raster()
  + facet_wrap(~test_prop,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_fill_viridis_c(expand=c(0,0))
  + ggtitle("log-space diff between simp and log_simp")
)
# ggsave("log-diff_simp-log.png",plot=fig_logdiff_simp_log, path = "./pix", width=3200,height=1800,units="px")

fig_logdiff_est_log <-(
  ggplot(dd,aes(phi,inc,fill=diff_est_log,group=test_prop))
  + geom_raster()
  + labs(x=bquote("Testing Focus Parameter"~phi), y=bquote("Prevalence"~Y), fill=bquote("log-Difference of \nExpected Positivity"~bar(P)))
  + facet_wrap(~test_prop,scale="free",labeller = label_bquote("Testing Proportion"~T~"="~.(test_prop)))
  # + scale_y_log10()
  + scale_fill_viridis_c(expand=c(0,0))
  + ggtitle("Numerical Underflow: log-space Diff between Linear Est and R qbeta output")
  + theme(axis.title.x = element_text(size = 18), # X-axis title font size
          axis.text.x = element_text(size = 16), # X-axis label font size
          axis.title.y = element_text(size = 18), # Y-axis title font size
          axis.text.y = element_text(size = 16), # Y-axis label font size
          plot.title = element_text(size = 18), # Plot title font size
          legend.title = element_text(size = 18),
          strip.text = element_text(size = 16))
)
fig_logdiff_est_log
ggsave("log-diff_est-log.png",plot=fig_logdiff_est_log, path = "./pix", width=4000,height=2000,units="px")

#### need to investigate more: what happens to a phi=0.99 case  
prop_pos_test_new(0.5,0.001,0.97,method = "cdf",debug = T)
prop_pos_test_new(0.5,0.001,0.97,method = "log",debug = T)
prop_pos_test_new(0.5,0.001,0.97,method = "simp",debug = T)

val_log <- function(i,phi,t) {
  phi_0 <- phi
  phi <- -log(1-phi)
  a <- i/phi; b <- (1-i)/phi
  qb <-qbeta(t,a,b, lower.tail = FALSE)
  log_qb <- log(qb)
  simp<-a/(a+b)*(t+(exp(a*log(qb)+b*log(1-qb)))/(beta(a,b)*a))
  log_simp <- log(simp)
  cdf <- pbeta(qb,a+1,b,lower.tail=FALSE)*(a/(a+b))
  log_cdf <- log(cdf)
  return(c(qb,log_qb,log_simp,log_cdf))
}

prop_pos_test_new(0.5,0.001,0.87,method = "log",debug = T)

val_log(0.5,0.93,0.001)
val_log(0.5,0.94,0.001)
val_log(0.5,0.95,0.001)

exp(val_log(0.8,0.97,0.001)[4])

####### Numerical observation
## even lwr output is "1", log(qbeta) could still have difference
## using both simp and log simp could somehow catch the difference
## However, for phi=0.99 very extreme case, even log(qbeta) gets solid zero

## very extreme case have positivity=io<since the log(1-qb)=-Inf
## cdf method will give "solid" zero for that case (log=-Inf)

dd_est <- (expand.grid(test_prop=c(0.05),#,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                   phi=seq(from=0.001,to=0.999,by=0.002),
                   inc=seq(from=0.001,to=0.999,by=0.002))
       %>% as_tibble()
       #%>% mutate(pos_prop_int=prop_pos_test_new(inc,test_prop,phi,method="int"))
       %>% mutate(pos_prop_est=prop_pos_test_new(inc,test_prop,phi,method="est"))
       %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
       #%>% mutate(diff_cdf_log=sign(pos_prop_cdf-pos_prop_log))
       %>% mutate(diff_est_log=abs_log_sub(pos_prop_est,pos_prop_log))
       )


fig_HR_logdiff_est_log<-(
  ggplot(dd_est,aes(phi,inc,fill=diff_est_log,group=test_prop))
  + geom_raster()
  + facet_wrap(~test_prop,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_fill_viridis_c(expand=c(0,0))
  + ggtitle("log-space diff between est and log_simp")
)
ggsave("NEW-log-diff_est-log.png",plot=fig_HR_logdiff_est_log, path = "./pix", width=3200,height=1800,units="px")

which(dd_est$inc==0.851)

slice_inc <- subset(dd_est, dd_est$inc==0.851)
slice_inc[which(is.nan(slice_inc$diff_est_log)),]
# print(slice_inc,n=1000)
# slice_inc 

dtF <- rbind(
  data.frame(phi=slice_inc$phi, y=slice_inc$pos_prop_log, class="logsimp_method", gp="pos_prop"),
  data.frame(phi=slice_inc$phi, y=slice_inc$pos_prop_est, class="est_method", gp = "pos_prop"),
  data.frame(phi=slice_inc$phi, y=slice_inc$diff_est_log, class="log_diff", gp="log_diff"))


fig_slice_logdiff_est_log <- (
  ggplot(data = dtF, mapping = aes(x = phi, y = y)) 
  +facet_grid(gp~., scale = "free_y")
  +geom_point(data=dtF[dtF$gp=="log_diff",], aes(color=class),size =0.6)
  +geom_point(data=dtF[dtF$gp=="pos_prop",], aes(color=class),size=0.4)
)
ggsave("slice_log_diff_est-logsimp.png",plot=fig_slice_logdiff_est_log, path = "./pix", width=3200,height=1800,units="px")


# val <- function(i,phi,t) {
#   phi_0 <- phi
#   phi <- -log(1-phi)
#   a <- i/phi; b <- (1-i)/phi
#   lwr <-qbeta(t,a,b, lower.tail = FALSE)
#   log_lwr <- log(lwr)
#   log_val = a/(a+b)*(t+(exp(a*log(lwr)+b*log(1-lwr)-log(beta(a,b)*a))))/t
#   est_val = (b+1-(1-lwr)*a*b)/(b+1-(1-lwr)*(a-1)*b)
#   est_num = (b+1-(1-lwr)*a*b)
#   est_denom = (b+1-(1-lwr)*(a-1)*b)
#   return(c(lwr,log_lwr,log_val,est_val,est_num,est_denom))
# }
# 
# ind<- 528
# val(as.numeric(slice_inc[ind,3]),as.numeric(slice_inc[ind,2]),as.numeric(slice_inc[ind,1]))
