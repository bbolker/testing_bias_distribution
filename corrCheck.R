library(purrr)
library(dplyr)

num_cor <- function(i, t, phi, n = 1e5) {
    gamma <- -log(1-phi)
    a <- i/gamma; b <- (1-i)/gamma
    lwr <- qbeta(t, a, b, lower.tail = FALSE)
    x <- rbeta(n, a, b)
    y <- as.numeric(x>lwr)
    cc <- cor.test(x,y)
    data.frame(gam=gamma, cor = cc$estimate, thresh=lwr)
}
phivec <- seq(0.1, 0.95, by = 0.05)
set.seed(101)
numcor <- (map_dfr(phivec, \(p) num_cor(i=0.1, t=0.01, phi=p))
    |> mutate(phi=phivec, .before = 1))

print(numcor)

