library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
source("testing_funs.R")

dd <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                   phi=seq(from=0.01,to=0.99,by=0.01),
                   inc=c(0.001,0.005,0.01,0.02,0.05,0.10,0.25,0.5))
       %>% as_tibble()
       %>% mutate(pos_prop=prop_pos_test_new(inc,test_prop,phi,method="log"))
)
# dd

print(ggplot(dd,aes(phi,pos_prop,col=log10(inc),group=test_prop))
      + geom_point(size=0.5)
      + facet_wrap(~test_prop,scale="free",labeller = label_both)
      #+ scale_y_log10()
      + scale_colour_viridis_c()
)
