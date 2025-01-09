library(dplyr)
library(ggplot2); theme_set(theme_classic(base_size=15))

library(shellpipes)
loadEnvironments()

base <- (ggplot(long)
	+ aes(date, count, color=category)
	+ geom_line()
	+ facet_wrap(~ pathogen, nrow=2)
	+ scale_color_brewer(palette="Dark2")
	+ guides(colour="none")
	+ theme(panel.grid.major = element_line())
)

summary(long)

## print(base + scale_y_log10() + facet_wrap(~ pathogen, scales="free"))
## print(base + facet_wrap(~ pathogen, scales="free"))

print(base + scale_y_log10())
## summary(dis)

(base + scale_y_log10()
	+ geom_line(data=dis, aes(date, percentP), color="blue")
) -> combplot

print(combplot %+% (long |> filter(category=="confirmations"))
	+ ylab("Confirmations / %Pos")
)

