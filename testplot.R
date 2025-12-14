library(shellpipes)
manageConflicts()

library(tidyr)
library(ggplot2); theme_set(theme_bw())
startGraphics()

wide <- rdsRead()
summary(wide)

long <- (wide
	|> pivot_longer(-t, names_to="group", values_to="prop")
)

summary(long)

print(ggplot(long)
	+ aes(t, prop, color = group)
	+ geom_line()
)
