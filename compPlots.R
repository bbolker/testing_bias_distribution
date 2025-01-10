library(shellpipes)
library(dplyr)
library(ggplot2); theme_set(theme_classic(base_size=15))

dat <- rdsRead()

testCut <- 0.1

summary(dat)

plist <- dat |> pull(prev) |> unique()

posPlot <- (ggplot(dat)
	+ aes(tests, positivity, shape=as.factor(add_lo), color=as.factor(prev))
	+ geom_line()
	+ scale_y_continuous(trans="logit")
	+ scale_y_continuous(trans="logit", breaks=c(0.01, 0.03, 0.1, 0.3, 0.5))
	+ guides(colour="none")
	+ scale_color_brewer(palette="Dark2")
	+ geom_hline(yintercept=plist, lty=3)
)

casePlot <- (ggplot(dat)
	+ aes(tests, posTests, shape=as.factor(add_lo), color=as.factor(prev))
	+ geom_line()
	+ scale_y_continuous(trans="logit", breaks=c(0.005, 0.01, 0.02, 0.03))
	+ guides(colour="none")
	+ scale_color_brewer(palette="Dark2")
	+ geom_hline(yintercept=plist, lty=3)
)

print(posPlot)
print(posPlot + xlim(NA, testCut))
print(casePlot)
print(casePlot + xlim(NA, testCut))
