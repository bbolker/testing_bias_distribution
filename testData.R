library(splines)
library(glmmTMB)
library(scam)

set.seed(1001)
n <- 1000
a <- 1
b <- 3
c <- 1
d <- 2
x <- sort(rexp(n))
v <- a + b/x
m <- c/(1+d*x^2)

y <- rnorm(n, mean=m, sd=sqrt(v))
dd <- data.frame(x, y)

gtfit <- glmmTMB(y~ns(x, 3)
               , dispformula = ~ns(x, 3)
               , data = dd
)
summary(gtfit)

gtfit2 <- update(gtfit,
                 dispformula = ~s(x, bs = "tp", k = 10),
                 REML = TRUE)

ypp <- predict(gtfit, se.fit=TRUE)

plot(x, ypp$fit, type="l")
lines(x, ypp$fit-ypp$se.fit, lty=2)
lines(x, ypp$fit+ypp$se.fit, lty=2)
points(x, y)

pred_se <- predict(gtfit, type = "disp")
pred_se2 <- predict(gtfit2, type = "disp")
plot(x, pred_se, type = "l", lwd =2)
lines(x, sqrt(v), col = 2, lwd = 2)
points(x, abs(residuals(gtfit)),
       col = adjustcolor("black",  alpha.f = 0.2))
lines(x, pred_se2, col = 4, lwd = 2)
