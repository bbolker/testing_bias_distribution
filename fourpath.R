library(shellpipes)
manageConflicts()

library(dplyr)
library(tidyr)

dat <- csvRead()

dat |> mutate_if(is.character, as.factor) |> summary()

nat <- (dat
	|> filter(geo_type=="nation")
	|> select(date=time_value, 
   	, sarscov2_tests, sarscov2_positive_tests
   	, hcov_tests, hcov_positive_tests
   	, flu_tests, flu_positive_tests=flua_positive_tests
   	, rsv_tests, rsv_positive_tests
	)
)

## Stretch to separate
long <- (nat
	|> pivot_longer(-date
		, values_to="count"
	) |> mutate(
		name = sub("positive_tests", "positives", name)
	) |> separate(name, "_", into = c("pathogen", "category"))
)
summary(long |> mutate_if(is.character, as.factor))
rdsSave(long)

## Now collapse to get positivity!
dis <- (long
	|> pivot_wider(names_from=category, values_from=count)
	|> mutate(V = positives/tests)
)
summary(dis |> mutate_if(is.character, as.factor))

rdaSave(long, dis)

