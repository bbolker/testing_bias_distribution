library(shellpipes)

library(dplyr)

base_lo <- seq(-10, 6, length.out=201)
prev <- 0.02
add_lo <- 5

bind_rows(NULL
	, tibble(base_lo, prev=0.01, add_lo=9)
	, tibble(base_lo, prev=0.0155, add_lo=5)
	, tibble(base_lo, prev=0.05, add_lo=3)
) |> mutate(NULL
	, posTests = prev*plogis(base_lo+add_lo)
	, negTests = (1-prev)*plogis(base_lo)
	, tests = posTests + negTests
	, positivity = posTests/tests
) |> rdsSave()

