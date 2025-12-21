library(shellpipes)
manageConflicts()

loadEnvironments()

summary(epi)
summary(obs)

library(ggplot2); theme_set(theme_bw())
startGraphics()

print(ggplot(epi)
	+ geom_line(aes(t, I/N))
	+ geom_line(aes(t, base), color="green")
	+ ggtitle("Prevalence and concern")
)

print(ggplot(obs)
	+ geom_line(aes(t, pos), color="red")
	+ geom_line(aes(t, neg), color="blue")
	+ ggtitle("tests")
)
