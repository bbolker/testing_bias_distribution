library(shellpipes)
manageConflicts(add="dplyr")

library(dplyr)
library(bbmle)

library(ggplot2); theme_set(theme_bw(base_size=15))
startGraphics()

loadEnvironments()

epi <- simulate(sirRates 
	, states = (list(t=0, S=N-I0, I=I0))
	, params = (list(beta=beta, D=D, N=N, deltat = deltat))
	, steps=steps
) |> mutate(series=1)


## We want to use rho=1 _unless_ it's being estimated
## since floating fits don't actually use rho
rho <-1
## Swap in the fitted parameters for the true ones
list2env(as.list(coef(rdsRead())), envir=.GlobalEnv)

new <- simulate(sirRates 
	, states = (list(t=0, S=N-I0, I=I0))
	, params = (list(beta=beta, D=D, N=N, deltat = deltat))
	, steps=steps
) |> mutate(series=2, I=I/rho)

comp <- (bind_rows(epi, new) 
	|> mutate(series=factor(series, labels=c("true", "inferred")))
)

eplot <- (ggplot(comp)
	+ aes(t, I, color=series)
	+ geom_line()
	+ theme(legend.position = c(0.02, 0.98)
		, legend.justification = c("left", "top")
	)
	+ ylab("Active cases")
)

print(eplot + comp |> filter(series=="true"))
print(eplot)
