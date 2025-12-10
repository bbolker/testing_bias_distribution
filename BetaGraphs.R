## This is a base-plot, beta function
bpic <- function(Y, T, phi, xtext=0.9, ytext=0.5, yoff=0.1){
  P = Y
  kap = log(1-phi)/(P*(P-1))
  
  alp = 1/(kap*(1-P))
  bet = 1/(kap*P)
  
  lwr <- qbeta(1-T, alp, bet)
  posTests <- integrate(\(x) dbeta(x, alp, bet)*x, lwr, 1)$value
  V <- posTests/T
  
  par(las=1, yaxs = "i",xaxs = "i")
  cc <- curve(  dbeta(x, alp, bet), from =0
              , to = 1
              , xlab = "x: probability of infection"
              , ylab = bquote(f[beta](x)~"prob density")
              , main = bquote("Y"==.(Y)~varphi==.(phi))
              )
  cc2 <- curve(dbeta(x, alp, bet), from = lwr, to = 1, add = TRUE)
  inf_ind <- which(cc2$y==Inf)
  cc2$y[inf_ind] <- 0
  polygon(c(cc2$x, rev(cc2$x)), c(rep(0, length(cc2$x)), rev(cc2$y)), col = "gray")
  
  height <- dbeta(P, alp, bet)
  
  text(x=xtext, y=(ytext+(1:-1)*yoff)*height , pos=2, labels=c(
      paste0("T = ", sprintf("%4.3f", T))
    , paste0("P = ", sprintf("%4.3f", V))
    , ""
  ))
  
  abline(v=P, lty=2)
}

bpic(Y=0.25, T=0.15, phi=0.01, yoff=0.15)
bpic(Y=0.25, T=0.15, phi=0.95, xtext=0.98, ytext=3.5, yoff=2)



library(ggplot2)
library(dplyr)
dd <- (expand.grid(par=seq(from=-4,to=4,by=0.01)
                   )
        %>% as_tibble()
        %>% mutate(hazard=1-exp(-par))
        %>% mutate(loghazard=1-exp(-exp(par)))
        %>% mutate(logOdd=exp(par)/(1+exp(par)))
        #%>% mutate(ratio=test_prop/pos_prop)
)
#dd

fig_hazard <- (
  ggplot(dd)
  + theme_bw()
  + geom_line(aes(par,hazard,col="Hazard"))
  + geom_line(aes(par,loghazard,col="logHazard"))
  #+ geom_line(aes(par,logOdd,col="logOdds"))
  + labs(x="Parameter", y="Probability", col="Approach")
  + ylim(0,1)
  + theme(axis.title.x = element_text(size = 18), # X-axis title font size
          axis.text.x = element_text(size = 14), # X-axis label font size
          axis.title.y = element_text(size = 18), # Y-axis title font size
          axis.text.y = element_text(size = 14), # Y-axis label font size
          #plot.title = element_text(size = 18), # Plot title font size
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 14)
          )
  #+ facet_wrap(~phi,scale="free",labeller = label_bquote(phi~"="~.(phi)))
  # + scale_y_log10()
  #+ scale_colour_viridis_c(trans="log10", breaks = brkvec)
  #+ ggtitle(bquote(bar(P)~" vs "~ Y ~", grouped by"~phi~", colored by"~T))
)

fig_hazard

Y <- 0.25
FocusTest <- function(x){
  if (x<=Y) {
    return(1)
  } else {
    return(Y/x)
  }
}
FocusTest <- Vectorize(FocusTest, vectorize.args = c("x"))

df <- (expand.grid(  Test=seq(from=0,to=1,by=0.01)
                   )
       %>% as_tibble()
       %>% mutate(Random = rep(Y, 101))
       %>% mutate(Focus=FocusTest(Test)))
df
plot(  df$Test
     , df$Random, type="l"
     , xlab = "Tested Proportion T"
     , ylab = bquote("Expected Positivity"~bar(P))
     #, main = "Totally Random Testing"
     )
text(x=0.2, y=Y, pos=3, labels=c(
  paste0("Y = ", sprintf("%4.2f", Y)))
  )
df$Focus

plot(  df$Test
     , df$Focus, type="l"
     , xlab = "Tested Proportion T"
     , ylab = bquote("Expected Positivity"~bar(P))
     #, main = "Perfectly Focused Testing"
     , ylim = c(0,1.2)
)
abline(v=Y, lty=2)
text(x=Y+0.2, y=0.1, pos=3, labels=c(
  paste0("Y = ", sprintf("%4.2f", Y)))
)

