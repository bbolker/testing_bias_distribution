## This is a base-plot, beta function
bpic <- function(Y, T, phi, xtext=0.9, ytext=0.5, yoff=0.1){
  P = Y
  kap = log(1-phi)/(P*(P-1))
  
  alp = 1/(kap*(1-P))
  bet = 1/(kap*P)
  
  lwr <- qbeta(1-T, alp, bet)
  posTests <- integrate(\(x) dbeta(x, alp, bet)*x, lwr, 1)$value
  V <- posTests/T
  
  #par(las=1, yaxs = "i",xaxs = "i")
  cc <- curve(  dbeta(x, alp, bet), from =0
              , to = 1
              , xlab = bquote(x~": probability of being infectious")
              , ylab = bquote(f[beta](x)~"probability density")
              #, main = bquote("Y"==.(Y)~phi==.(phi))
              )
  cc2 <- curve(dbeta(x, alp, bet), from = lwr, to = 1, add = TRUE)
  inf_ind <- which(cc2$y==Inf)
  cc2$y[inf_ind] <- 0
  polygon(c(cc2$x, rev(cc2$x)), c(rep(0, length(cc2$x)), rev(cc2$y)), col = "gray")
  
  height <- dbeta(P, alp, bet)
  
  text(x=xtext, y=(ytext+(1:-1)*yoff)*height , pos=2, labels=c(
      paste0("T = ", sprintf("%4.3f", T))
    , paste0("E[P] = ", sprintf("%4.3f", V))
    , ""
  ))
  
  abline(v=P, lty=2)
}
# ggsave("test_positivity_vs_phi-test_prop.png",plot=fig_pos_vs_phi, path = "./pix", width=3200,height=1800,units="px")

bpic(Y=0.25, T=0.1, phi=0.01, yoff=0.15)
bpic(Y=0.25, T=0.15, phi=0.95, xtext=0.98, ytext=3.5, yoff=2)




# bpic(P=0.25, T=0.1, kap=0.1)
# bpic(P=0.25, T=0.2, kap=0.1)
# bpic(P=0.25, T=0.1, kap=1)
# bpic(P=0.25, T=0.2, kap=1)

# plot(c(1, 9), 1:2, type = "n")
# polygon(1:9, c(2,1,2,1,NA,2,1,2,1),
#         col = c("red", "blue"),
#         border = c("green", "yellow"),
#         lwd = 3, lty = c("dashed", "solid"))
# ?rev




