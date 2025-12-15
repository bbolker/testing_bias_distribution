library(shellpipes)
manageConflicts()

library(ggplot2); theme_set(theme_bw())
startGraphics()

print(ggplot(rdsRead())
	+ aes(t, pos)
	+ geom_line()
)
