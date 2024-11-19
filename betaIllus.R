library(shellpipes)

bpic <- function(P, T, kap, xtext=0.9, ytext=0.5, yoff=0.1){

	alp = 1/(kap*(1-P))
	bet = 1/(kap*P)

	lwr <- qbeta(1-T, alp, bet)
	posTests <- integrate(\(x) dbeta(x, alp, bet)*x, lwr, 1)$value
	V <- posTests/T

	par(las=1, yaxs = "i",xaxs = "i")
	cc <- curve(dbeta(x, alp, bet), from =0, to = 1, xlab = "prob(infected)", ylab = "prob density")
	cc2 <- curve(dbeta(x, alp, bet), from = lwr, to = 1, add = TRUE)
	polygon(c(cc2$x, rev(cc2$x)), c(rep(0, length(cc2$x)), rev(cc2$y)), col = "gray")

	height <- dbeta(P, alp, bet)

	text(x=xtext, y=(ytext+(1:-1)*yoff)*height , pos=2, labels=c(
		paste0("T = ", sprintf("%4.3f", T))
		, paste0("V = ", sprintf("%4.3f", V))
		, paste0("C = ", sprintf("%4.3f", T*V))
	))
}

bpic(P=0.25, T=0.1, kap=0.1)
bpic(P=0.25, T=0.2, kap=0.1)
bpic(P=0.25, T=0.1, kap=1)
bpic(P=0.25, T=0.2, kap=1)
