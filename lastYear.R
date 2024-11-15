library(shellpipes)
manageConflicts()

library(dplyr)

dat <- csvRead()

## dat |> mutate_if(is.character, as.factor) |> summary()

nat <- (dat
	|> filter(geo_type=="nation")
	|> select(time_value, 
   	, sarscov2_tests  sarscov2_positive_tests
   	, hcov_tests      hcov_positive_tests flu_positive_tests
   	, flu_tests     fluah1pdm09_positive_tests fluah3_positive_tests
	)
) 

