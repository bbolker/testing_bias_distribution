# Sys.setenv(LANG = "en")

library(shellpipes)
rpcall("const.fixed.check.Rout check.fit.R const.fixed.fit.rda")
rpcall("const.flex.check.Rout check.fit.R const.flex.fit.rda")
rpcall("float.flex.check.Rout check.fit.R float.flex.fit.rda")

suppressPackageStartupMessages(library(dplyr))

suppressPackageStartupMessages(library(macpan2))
options(macpan2_verbose = FALSE)

library(tidyr)
library(ggplot2); theme_set(theme_bw())
loadEnvironments()

(sir_sim
  |> mp_tmb_update(
    default = fit_bklist
  )
  |> mp_simulator(
    time_steps = 30
    , outputs = c("OT","OP","I","S","B_lik")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  |> pivot_wider(names_from = matrix,values_from = value)
) -> check_fit

loglik<-fit$objective

startGraphics()
fit_curve <- (ggplot() + theme_bw()
+ geom_point(data = dat_pall, aes(time,I,color="I(t)",shape="Data"))
+ geom_point(data = dat_pall, aes(time,OT,color="OT(t)",shape="Data"))
+ geom_point(data = dat_pall, aes(time,OP,color="OP(t)",shape="Data"))
+ geom_line(data = check_fit, aes(time+tmin-1,I,color="I(t)", linetype ="Fitted Model Simulation"))
+ geom_line(data = check_fit, aes(time+tmin-1,OT,color="OT(t)", linetype ="Fitted Model Simulation"))
+ geom_line(data = check_fit, aes(time+tmin-1,OP,color="OP(t)", linetype ="Fitted Model Simulation"))
+ labs(x="Time t", y="Case Count")
+ ylim(0,35000)
+ labs(title = bquote("loglik="~.(loglik)))
)
print(fit_curve)
#### Not converge (over iterate/function limit) for small population (1e-5) and limited data (t=10-20)
#### identifiability issue with large population and larger data set.
