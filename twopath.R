library(shellpipes)
manageConflicts()

library(dplyr)
library(tidyr)

dat <- csvRead()

dat |> mutate_if(is.character, as.factor) |> summary()

nat <- (dat
	|> filter(geo_type=="nation")
	|> select(date=time_value, 
   	, flu_tests, flu_positive_tests=flua_positive_tests
   	, rsv_tests, rsv_positive_tests
	)
)

## Stretch to separate
long <- (nat
	|> pivot_longer(-date
		, values_to="count"
	) |> mutate(
		name = sub("positive_tests", "confirmations", name)
	) |> separate(name, "_", into = c("pathogen", "category"))
)
summary(long |> mutate_if(is.character, as.factor))
rdsSave(long)

## Now collapse to get positivity!
dis <- (long
	|> pivot_wider(names_from=category, values_from=count)
	|> mutate(percentP = 100*confirmations/tests)
)
summary(dis |> mutate_if(is.character, as.factor))

## And re-expand?
new <- (dis
	|> pivot_longer(-c(date, pathogen)
		, values_to="count"
	)
)
summary(new |> mutate_if(is.character, as.factor))
rdaSave(long, dis, new)

