library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
source("testing_funs.R")

dd <- (expand.grid(test_prop=c(0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.25,0.5),
                   Phi=seq(from=0.001,to=0.999,by=0.001),
                   prev=c(0.001,0.005,0.01,0.02,0.05,0.10,0.25,0.5))
       %>% as_tibble()
       %>% mutate(B=Bfun_HR(prev,test_prop,Phi))
       %>% mutate(T_B=1-B)
       %>% mutate(T_V=1-B*Phi)
       %>% mutate(pos_prop=Pfun_HR(prev,test_prop,Phi))
)
# dd[which(dd$B>1),]
# which.max(dd$B)
# dd[6238,]
# as.numeric(dd[946,])
# (1-0.001)/(1-0.005*(1-0.07))

brkvec <- c(10^(-3:-1), 0.5)

figHR_pos_vs_phi <- (
  ggplot(dd,aes(Phi,pos_prop,col=prev,group=test_prop))
  + geom_point(size=0.5)
  + facet_wrap(~test_prop,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("test positivity vs phi, grouped by testing proportion, colored by prevalence")
)
ggsave("HR_test_positivity_vs_phi-test_prop.png",plot=figHR_pos_vs_phi, path = "./pix", width=3200,height=1800,units="px")

# figHR_TB_vs_phi <- (
#   ggplot(dd,aes(Phi,T_B,col=prev,group=test_prop))
#   + geom_point(size=0.5)
#   + facet_wrap(~test_prop,scale="free",labeller = label_both)
#   # + scale_y_log10()
#   + scale_colour_viridis_c(trans="log10", breaks = brkvec)
#   + ggtitle("T_B vs phi, grouped by testing proportion, colored by prevalence")
# )

# figHR_TB_vs_phi

# figHR_TV_vs_phi <- (
#   ggplot(dd,aes(Phi,T_V,col=prev,group=test_prop))
#   + geom_point(size=0.5)
#   + facet_wrap(~test_prop,scale="free",labeller = label_both)
#   # + scale_y_log10()
#   + scale_colour_viridis_c(trans="log10", breaks = brkvec)
#   + ggtitle("T_V vs phi, grouped by testing proportion, colored by prevalence")
# )
# 
# figHR_TV_vs_phi

dd2 <- (expand.grid(Phi=c(0.01,0.02,0.05,0.1,0.2,0.5,0.7,0.9,0.99),
                    test_prop=seq(from=0.001,to=0.999,by=0.001),
                    prev=c(0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.20,0.5))
        %>% as_tibble()
        %>% mutate(B=Bfun_HR(prev,test_prop,Phi))
        %>% mutate(T_B=1-B)
        %>% mutate(T_V=1-B*Phi)
        %>% mutate(pos_prop=Pfun_HR(prev,test_prop,Phi))
        )

figHR_pos_vs_test_prop <- (
  ggplot(dd2,aes(test_prop,pos_prop,col=prev,group=Phi))
  + geom_point(size=0.5)
  + facet_wrap(~Phi,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("test Positivity vs test proportion, grouped by phi, colored by prevalence")
)
ggsave("HR_test_positivity_vs_test_proportion-phi.png",plot=figHR_pos_vs_test_prop, path = "./pix", width=3200,height=1800,units="px")


#### Heatmap
dd3 <- (expand.grid(Phi=seq(from=0.001,to=0.999,by=0.001),
                    test_prop=seq(from=0.001,to=0.999,by=0.001),
                    prev=c(0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.20,0.5))
        %>% as_tibble()
        %>% mutate(pos_prop=Pfun_HR(prev,test_prop,Phi))
)

figHR_heatmap <- (
  ggplot(dd3,aes(test_prop,Phi,fill=pos_prop,group=prev))
  + geom_raster()
  + facet_wrap(~prev,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_fill_viridis_c()
  + ggtitle("HR Heatmap: test proportion vs Phi, grouped by prev, colored by positivity")
)
ggsave("HR_heatmap.png",plot=figHR_heatmap, path = "./pix", width=3200,height=1800,units="px")
