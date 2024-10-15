library(tidyverse)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(bbmle)
library(DPQ)
library(deSolve)
source("testing_funs.R")

########## incident/testing positivity ratio
## Temporarily remove overflow case caused by lwr=1 for better ilustration
HL_fn <- function(prop_pos,i,t,phi){
  phi<- -log(1-phi)
  a <- i/phi; 
  b <- (1-i)/phi
  lwr <- qbeta(t,a,b,lower.tail=FALSE)
  if (log(lwr)==0){
    return(NaN)
  }
  return(i/prop_pos)
}
HL_fn <- Vectorize(HL_fn,c("prop_pos","i","t","phi"))

dd_ratio <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                         phi=seq(from=0.01,to=0.99,by=0.01),
                         inc=seq(from=0.01,to=0.99,by=0.01))
             %>% as_tibble()
             %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
             %>% mutate(ratio=HL_fn(pos_prop_log,inc,test_prop,phi))
)

print(ggplot(dd_ratio,aes(phi,inc,fill=ratio,group=test_prop))
      + geom_tile()
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_fill_viridis_c()
      + ggtitle("inc/pos_prop ratio, group by test_prop")
)

##### Focus on more reasonable intervals for incidence
dd_ratio_focus <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                               phi=seq(from=0.001,to=0.999,by=0.001),
                               inc=seq(from=0.001,to=0.25,by=0.001))
                   %>% as_tibble()
                   %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
                   %>% mutate(ratio=HL_fn(pos_prop_log,inc,test_prop,phi))
)

print(ggplot(dd_ratio_focus,aes(phi,inc,fill=log(ratio),group=test_prop))
      + geom_tile()
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_fill_viridis_c()
      + ggtitle("inc/pos_prop ratio, group by test_prop")
)

##### Slice by different inc
inc_slice<-(expand.grid(test_prop=c(0.001,0.002,0.005,0.010,0.020,0.050,0.1,0.2,0.5),
                        phi=seq(from=0.001,to=0.999,by=0.001),
                        inc=c(0.001,0.002,0.005,0.010,0.020,0.050,0.1,0.2,0.5))
            %>% as_tibble()
            %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
            %>% mutate(ratio=HL_fn(pos_prop_log,inc,test_prop,phi))
)

print(ggplot(inc_slice,aes(phi,ratio,color=log(test_prop),group=inc))
      + geom_point(size=0.5)
      + facet_wrap(~inc,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_colour_viridis_c()
      + ggtitle("inc/pos_prop ratio, slice at different incidence")
)

###### Investigation: what happens to the decrease on yellow curves?
print(inc_slice[which(inc_slice$test_prop==.5 & inc_slice$inc==0.001),],n=1000)
### pos_prob start to increase suddenly

t1 <- inc_slice[which(inc_slice$test_prop==.5 & inc_slice$inc==0.001),]

ii1 <- as.numeric(t1[500,3])
phi1 <- as.numeric(t1[500,2])
aa1 <- ii1/phi1
bb1 <- (1-ii1)/phi1


ii2 <- as.numeric(t1[900,3])
phi2 <- as.numeric(t1[900,2])
aa2 <- ii2/phi2
bb2 <- (1-ii2)/phi2

plot(seq(0.001,0.999,0.001),dbeta(seq(0.001,0.999,0.001),aa1,bb1),type="l")
points(seq(0.001,0.999,0.001),dbeta(seq(0.001,0.999,0.001),aa2,bb2),type="l",col="red")

###### Investigation: what happens to the decrease on large test_prop with large inc?
print(inc_slice[which(inc_slice$test_prop==0.001 & inc_slice$inc==0.5),],n=1000)
### Overflow
### FIXME???
t2 <- inc_slice[which(inc_slice$test_prop==0.001 & inc_slice$inc==0.5),]

prop_pos_test_new(t2[934,3], t2[934,1], t2[934,2], debug=T)

prop_pos_test_new(t2[935,3], t2[935,1], t2[935,2], debug=T)
prop_pos_test_new(t2[935,3], t2[935,1], t2[935,2], method="cdf",debug=T)

prop_pos_test_new(t2[936,3], t2[936,1], t2[936,2], debug=T)
prop_pos_test_new(t2[936,3], t2[936,1], t2[936,2], method="cdf",debug=T)

prop_pos_test_new(t2[940,3], t2[940,1], t2[940,2], debug=T)
prop_pos_test_new(t2[940,3], t2[940,1], t2[940,2], method="cdf",debug=T)

## cdf method even worse
## A "good" example of log-simp is slightly better?

##### Slice by different phi
phi_slice <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                          phi=seq(from=0.1,to=0.9,by=0.1),
                          inc=seq(from=0.001,to=0.50,by=0.001))
              %>% as_tibble()
              %>% mutate(pos_prop_log=prop_pos_test_new(inc,test_prop,phi,method="log"))
              %>% mutate(ratio=HL_fn(pos_prop_log,inc,test_prop,phi))
)

print(ggplot(phi_slice,aes(inc,ratio,color=phi,group=test_prop))
      + geom_point(size=0.5)
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_colour_viridis_c()
      + ylim(0, 0.8)
      + ggtitle("inc/pos_prop ratio, slice at different phi and test_prop")
)
