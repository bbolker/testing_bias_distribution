# Sys.setenv(LANG = "en")

library(shellpipes)
rpcall("const.flex.pois.check.Rout check.pois.fit.R const.flex.pois.fit.rda")
rpcall("const.fixed.pois.check.Rout check.pois.fit.R const.fixed.pois.fit.rda")
rpcall("float.flex.pois.check.Rout check.pois.fit.R float.flex.pois.fit.rda")

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
    , outputs = c("OPos","ONeg","I","S","B_lik")
  ) 
  |> mp_trajectory()
  |> dplyr::select(-c(row, col))
  |> pivot_wider(names_from = matrix,values_from = value)
) -> check_fit

loglik<-fit$objective

startGraphics()
fit_curve <- (ggplot() + theme_bw()
+ geom_point(data = dat_pall, aes(time,I,color="I(t)",shape="Data"))
+ geom_point(data = dat_pall, aes(time,OPos,color="OPos(t)",shape="Data"))
+ geom_point(data = dat_pall, aes(time,ONeg,color="ONeg(t)",shape="Data"))
+ geom_line(data = check_fit, aes(time+tmin-1,I,color="I(t)", linetype ="Fitted Model"))
+ geom_line(data = check_fit, aes(time+tmin-1,OPos,color="OPos(t)", linetype ="Fitted Model"))
+ geom_line(data = check_fit, aes(time+tmin-1,ONeg,color="ONeg(t)", linetype ="Fitted Model"))
+ labs(x="Time t", y="Case Count")
# + ylim(0,50000)
+ labs(title = bquote("-loglik="~.(-loglik)))
)
print(fit_curve)
#### Not converge (over iterate/function limit) for small population (1e-5) and limited data (t=10-20)
#### identifiability issue with large population and larger data set.
