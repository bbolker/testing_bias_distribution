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

rdsSave(grid)
