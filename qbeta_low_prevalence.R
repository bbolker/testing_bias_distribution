library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
source("testing_funs.R")

log(qbeta(0.2,1,2,lower.tail = FALSE)+qbeta(0.2,2,1))


prop_pos_test_new(0.001,0.4,0.99,method="est",debug = T)
prop_pos_test_new(0.001,0.572,0.99,method="est",debug = T)


Phi <- 0.99
i <- 0.001
t <- 0.57
phi <- -log(1-Phi)
a <- i/phi
b <- (1-i)/phi

plot(seq(0,1,0.001),dbeta(seq(0,1,0.001),a,b),type="l")


lwr <- qbeta(t,a,b,lower.tail=FALSE)
lwr
pbeta(lwr,a,b,lower.tail = FALSE)

g <- function(l,a,b){
  l^a/(a*(a+1))*(a+1-l*(b-1)*(a))
}
g(lwr,a,b)
g(lwr,a+1,b)

(beta(a+1,b)-g(lwr,a+1,b))/(beta(a,b)*t)
prop_pos_test_new(i,t,Phi,method="est",debug = T)


dd <- (expand.grid(test_prop=seq(from=0.001,to=0.999,by=0.001),
                   Phi=seq(from=0.001,to=0.999,by=0.001),
                   prev=c(0.001))
       %>% as_tibble()
       %>% mutate(phi = -log(1-Phi))
       %>% mutate(a= prev/phi)
       %>% mutate(b= (1-prev)/phi)
       %>% mutate(l=qbeta(test_prop,a,b,lower.tail=FALSE))
       )

fig_lwr_heatmap <- (
  ggplot(dd,aes(test_prop,Phi,fill=log(l)))
  + geom_raster()
  + scale_fill_viridis_c()
  + ggtitle("qbeta for prev=0.001 as a function of phi and test proportion")
)

ggsave("qbeta_heatmap_low_prevalence.png",plot=fig_lwr_heatmap, path = "./pix", width=3200,height=1800,units="px")


