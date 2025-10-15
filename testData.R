library(splines)
library(glmmTMB)

set.seed(1001)
n <- 100
a <- 1
b <- 1
c <- 1
d <- 2
x <- sort(rexp(n))
v <- a + b/x
m <- c/(1+d*x^2)

y <- rnorm(n, mean=m, sd=sqrt(v))

gtfit <- glmmTMB(y~ns(x, 3)
	, dispformula = ~ns(x, 3)
)
summary(gtfit)

ypp <- predict(gtfit, se.fit=TRUE)

plot(x, ypp$fit, type="l")
lines(x, ypp$fit-ypp$se.fit, lty=2)
lines(x, ypp$fit+ypp$se.fit, lty=2)
points(x, y)

