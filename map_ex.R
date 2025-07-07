## mapping/profiling
library(macpan2)
library(fitode)

sir <- (mp_tmb_library("starter_models", "sir", package = "macpan2")
  |> mp_tmb_insert(default = list(beta=70/52, gamma=60/52,
                                  N=40000, I=1, phi=6))
)

## convert data to macpan2 format
sirdat <- fitode::SierraLeone2014 |>
  dplyr::transmute(
    ## can't handle year + decimal time
    time = seq_along(times),
    matrix = "I",
    row = 0,
    col = 0,
    value = confirmed)

sir_calibrator = mp_tmb_calibrator(sir
  , data = sirdat,
  , traj = list(I = mp_neg_bin(disp=mp_fit(2))),
  , par = c("beta", "gamma")
)

mp_optimize(sir_calibrator)


proffun <- function(..., plot.it = TRUE) {
  pp <- mp_tmb_profile(sir_calibrator, "beta", trace = FALSE, ...)
  minval <- min(pp$value, na.rm = TRUE)
  pp <- transform(pp, value = value - minval)
  attr(pp, "minval") <- minval
  if (plot.it) plot(value ~ params, pp, type = "b")
  invisible(pp)
}
tmbobj <- mp_tmb(sir_calibrator)
mp_parameterization(sir_calibrator)
par(las=1, bty = "l")
proffun()

## can plot over a wider range:
pp4 <- proffun(ytol = 4)
## don't know why profile gets wonky there, or why upper end of params range doesn't go farther

pp5 <- proffun(ytol = Inf, parm.range = c(0.5, 4.5), maxit = 1e5)



