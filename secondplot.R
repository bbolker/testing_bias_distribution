library(dplyr)
library(ggplot2); theme_set(theme_classic(base_size=15))

library(shellpipes)
loadEnvironments()

summary(new)

base <- (ggplot(new |> filter(name != "tests"))
	+ aes(date, count, color=pathogen)
	+ geom_line()
	+ facet_wrap(~ name, nrow=2, scale="free")
	+ scale_color_brewer(palette="Dark2")
	+ ylab("")
	## + guides(colour="none")
)

print(base)

