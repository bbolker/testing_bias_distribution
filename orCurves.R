library(shellpipes)

library(dplyr)

base_lo <- seq(-6, 6, length.out=201)
prev <- c(0.01, 0.03)
add_lo <- 2*(1:3)

grid <- (expand.grid(base_lo=base_lo, prev=prev, add_lo=add_lo)
	|> mutate(NULL
		, posTests = prev*plogis(base_lo+add_lo)
		, negTests = (1-prev)*plogis(base_lo)
		, tests = posTests + negTests
		, positivity = posTests/tests
	)
)

print(grid)

library(ggplot2); theme_set(theme_classic(base_size=15))
posPlot <- (ggplot(grid)
	+ aes(tests, positivity, shape=as.factor(add_lo), color=as.factor(prev))
	+ geom_line()
	+ scale_y_continuous(trans="logit")
	+ scale_y_continuous(trans="logit", breaks=c(0.01, 0.03, 0.1, 0.3))
	+ guides(colour="none")
)

casePlot <- (ggplot(grid)
	+ aes(tests, posTests, shape=as.factor(add_lo), color=as.factor(prev))
	+ geom_line()
	+ scale_y_continuous(trans="logit", breaks=c(0.005, 0.01, 0.02, 0.03))
	+ guides(colour="none")
)

print(posPlot)
print(posPlot + xlim(NA, 0.2))
print(casePlot)
print(casePlot + xlim(NA, 0.2))

