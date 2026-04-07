# Sys.setenv(LANG = "en")

library(shellpipes)
rpcall("float.plot.Rout plot.dataview.R float.data.rda")
rpcall("const.plot.Rout plot.dataview.R const.data.rda")
rpcall("float.pois.plot.Rout pois.dataview.R float.pois.data.rda")
rpcall("const.pois.plot.Rout pois.dataview.R const.pois.data.rda")

suppressPackageStartupMessages(library(dplyr))
library(tidyr)
library(ggplot2); theme_set(theme_bw())
loadEnvironments()

startGraphics()
model_curve<-(ggplot() + theme_bw()
  + geom_line(data = dat_pall, aes(time,I,color="I(t)"))
  + geom_line(data = dat_pall, aes(time,ONeg,color="Neg(t)"))
  + geom_line(data = dat_pall, aes(time,OPos,color="Pos(t)"))
  + geom_point(data = dat_fit, aes(time,ONeg,color="Neg(t)",shape="Fitted data"))
  + geom_point(data = dat_fit, aes(time,OPos,color="Pos(t)",shape="Fitted data"))
  + labs(x="Time t", y="Case Count")
)
print(model_curve)
#ggsave("FloatSIR_ModelCurve.png",plot=model_curve, path = "../pix", width=1800,height=900,units="px")

## Compare true baseline T_B and theoretical likelihood 1-B_lik
Phi <- exp(-h)
neg <- c(dat_pall$ONeg)
pos <- c(dat_pall$OPos)

B_lik <- 1/(2*N*Phi)*(((N-neg)*Phi+N-pos)-sqrt(((N-neg)*Phi+N-pos)^2-4*N*Phi*(N-pos-neg)))
# dat_pall$T_B

TB_Curve<-(ggplot() + theme_bw()
              + geom_line(data = dat_pall, aes(time,T_B,color="T_B"))
              + geom_line(data = dat_pall, aes(time,1-B_lik,color="likelihood B- Sol"))              + labs(x="Time t", y="Probability")
)
print(TB_Curve)