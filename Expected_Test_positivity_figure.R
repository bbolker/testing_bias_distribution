library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(viridis)
library(glue)
source("testing_funs.R")

dd <- (expand.grid(test_prop=c(0.001,0.0012,0.0015,0.002,0.005,0.01,0.05,0.1,0.25),
                   phi=seq(from=0.01,to=0.99,by=0.01),
                   prev=c(0.001,0.005,0.01,0.02,0.05,0.10,0.25,0.5))
       %>% as_tibble()
       %>% mutate(pos_prop=prop_pos_test_new(prev,test_prop,phi,method="est"))
)
# dd

brkvec <- c(10^(-3:-1), 0.5)

fig_pos_vs_phi <- (
  ggplot(dd,aes(phi,pos_prop,col=prev,group=test_prop))
  + geom_point(size=0.5)
  + facet_wrap(~test_prop,scale="free",labeller = label_both)
  # + scale_y_log10()
  + scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle("test positivity vs phi, grouped by testing proportion, colored by prevalence")
)
fig_pos_vs_phi
ggsave("test_positivity_vs_phi-test_prop.png",plot=fig_pos_vs_phi, path = "./pix", width=3200,height=1800,units="px")

dd2 <- (expand.grid(phi=c(0.001,0.01,0.05,0.1,0.2,0.5,0.9,0.95,0.99),
                    test_prop=seq(from=0.01,to=0.99,by=0.01),
                    prev=c(0.001,0.005,0.01,0.02,0.05,0.10,0.25,0.5))
        %>% as_tibble()
        %>% mutate(pos_prop=prop_pos_test_new(prev,test_prop,phi,method="est"))
)


fig_pos_vs_test_prop <- (
  ggplot(dd2,aes(test_prop,pos_prop,col=prev,group=phi))
  + geom_point(size=0.5)
  + labs(x=bquote(T), y=bquote(bar(P)), col=bquote(Y))
  + facet_wrap(~phi,scale="free",labeller = label_bquote(phi~"="~.(phi)))
  # + scale_y_log10()
  + scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle(bquote(bar(P)~" vs "~ T ~", grouped by"~phi~", colored by"~Y))
)
fig_pos_vs_test_prop
ggsave("test_positivity_vs_test_proportion-phi.png",plot=fig_pos_vs_test_prop, path = "./pix", width=3200,height=1800,units="px")


dd3 <- (expand.grid(phi=c(0.001,0.01,0.1,0.99),
                    prev=seq(from=0.001,to=0.999,by=0.001),
                    test_prop=c(0.001,0.005,0.01,0.05,0.10,0.2,0.5,0.9))
        %>% as_tibble()
        %>% mutate(pos_prop=prop_pos_test_new(prev,test_prop,phi,method="est"))
        #%>% mutate(ratio=test_prop/pos_prop)
)


fig_pos_vs_prev <- (
  ggplot(dd3,aes(prev,pos_prop,col=test_prop,group=phi))
  + geom_point(size=0.5)
  + labs(x=bquote("Prevalence"~Y), y=bquote("Expected Positivit"~bar(P)), col=bquote("Testing Proportion"~T))
  + facet_wrap(~phi,scale="free",labeller = label_bquote("Testing Focus"~phi~"="~.(phi)))
  # + scale_y_log10()
  + ylim(0,1)
  + scale_colour_viridis_c(trans="log10", breaks = brkvec)
  + ggtitle(bquote("Expected Positivity"~bar(P)~" vs Prevalence"~ Y ~", Grouped by Testing Focus"~phi~", Colored by Testing Proportion"~T))
  )
fig_pos_vs_prev
ggsave("test_positivity_vs_prev-phi.png",plot=fig_pos_vs_prev, path = "./pix", width=3200,height=1800,units="px")
