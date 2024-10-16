library(shellpipes)
library(fitode)
if (packageVersion("fitode") <= "0.1.1") {
    stop("please install the latest version of fitode via remotes::install_github('parksw3/fitode')")
}
if (!require("RTMBode")) {
    stop(" please install RTMBode via remotes::install_github('kaskr/RTMB/RTMBode')")
}
library(RTMBode)
library(RTMB)

SierraLeone2014b <- rbind(
    c(times=SierraLeone2014$times[1] -
          diff(SierraLeone2014$times)[1], confirmed=NA),
    SierraLeone2014
)
SIR_model <- odemodel(
    name="SIR (nbinom)",
    model=list(
        S ~ - beta * S * I/N,
        I ~ beta * S * I/N - gamma * I,
        R ~ gamma * I
    ),
    observation=list(
        confirmed ~ dnbinom(mu=R, size=phi)
    ),
    initial=list(
        S ~ N * (1 - i0),
        I ~ N * i0,
        R ~ 0
    ),
    diffnames="R",
    par=c("beta", "gamma", "N", "i0", "phi"),
    link=c(i0="logit")
)

SIR_start <- c(beta=70, gamma=60, N=40000, i0=0.0004, phi=6)

set.seed(101)
ss_SIR <- simulate(SIR_model,
                   parms=SIR_start, times=SierraLeone2014b$times)

plot(SierraLeone2014)
lines(ss_SIR$times, ss_SIR$R)

## tuning the proposal distribution is an important step for
## efficiently sampling the posterior distribution
## https://jellis18.github.io/post/2018-01-02-mcmc-part1/
proposal.vcov <- matrix(0, 5, 5)
diag(proposal.vcov) <- c(1e-4, 1e-4, 1e-4, 1e-8, 1e-4)

## using very vague priors
## maybe a useful document
## http://www.stat.columbia.edu/~gelman/research/published/p039-_o.pdf
SIR_fit <- fitodeMCMC(
    model=SIR_model,
    data=SierraLeone2014b,
    start=SIR_start,
    chains = 1,
    iter = 200,
    burnin = 100,
    thin = 1,
    proposal.vcov=proposal.vcov,
    prior = list(
        beta ~ dgamma(shape=2, rate=1/30),
        gamma ~ dgamma(shape=2, rate=1/30),
        N ~ dgamma(shape=2, rate=1/20000),
        i0 ~ dbeta(shape1=4, shape2=9996),
        phi ~ dgamma(shape=2, rate=1/3)
    )
)

## accessing posterior
## not looking great because
## (1) chain too short
## (2) need to tune in proposal distribution
plot(SIR_fit@mcmc[[1]][,1])

confint(SIR_fit)

## before fixing fitode method, get NaN error in quantile function
## "missing values and NaN's not allowed if 'na.rm' is FALSE"
plot(SIR_fit, level=0.95)

predict(SIR_fit, level=0.95)

## attempt to reimplement all of this in RTMBode
library(RTMB)
## regular gradient function as you would use with deSolve()
SIRmod <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        incidence <- beta * S * I/N
        recovery <- gamma * I
        return(list(c(S = - incidence,
                      I = incidence - recovery,
                      R = recovery)))
    })
}

## link function
SIR_logstart <- log(SIR_start) |>
    setNames(paste0("log_", names(SIR_start)))

likfun <- function(pars) {
    ## inverse-link function
    for (nm in names(pars)) {
        assign(gsub("log_", "", nm), exp(pars[[nm]]))
    }
    ini <- c(S=N*(1-i0), I = N*i0, R = 0)
    ## get subset of *dynamical* parameters
    ## unlist(tibble::lst(beta, gamma, N)) (shortcut for self-named list) causes trouble
    ## or? use helper function  vdiff <- function(x, nm) {x[!names(x) %in% nm]}
    ode_pars <- c(beta = beta, gamma = gamma, N = N)
    sol <- ode(func = SIRmod,
               y = ini,
               parms = ode_pars,
               times = SierraLeone2014b$times)
    mu <- diff(sol[,"R"])
    ## base-R dnbinom() is difficult; use dnbinom2, specify var rather than phi ('size')
    var <- mu*(1+mu/phi)
    nll <- -sum(dnbinom2(SierraLeone2014b$confirmed[-1],
                         mu = mu , var = var, log = TRUE))
    ## negative log(prior): must use dgamma() with scale, not rate parameter
    nlog_prior <- -1*(
        dgamma(beta, shape=2, scale=30, log = TRUE) +
        dgamma(gamma, shape=2, scale=30, log = TRUE) +
        dgamma(N, shape=2, scale=20000, log = TRUE) +
        dbeta(i0, shape1=4, shape2=9996, log = TRUE) +
        dgamma(phi, shape=2, scale=3, log = TRUE))
    nll + nlog_prior
}

likfun(SIR_logstart)
ff <- MakeADFun(likfun, as.list(SIR_logstart))
ff$fn()
with(ff, nlminb(par, fn, gr))

## RTMB version is much faster, even without consider gradient-computation advantage,
## (which is irrelevant for Metropolis-Hastings)
microbenchmark::microbenchmark(
                    raw = likfun(SIR_logstart),
                    rtmb = ff$fn(ff$par))

## BUT, this doesn't work well ...

library(MCMCpack)
llfun <- function(x) -ff$fn(x)
res <- MCMCmetrop1R(llfun, theta.init=ff$env$last.par.best)
matplot(res, type = "l")

library(tmbstan)
## starting from 'random', or not specifying priors, breaks; very slow sampling
## eventually fails
## set bounds to try to prevent badness ...
lpb <- ff$env$last.par.best
tt <- tmbstan(ff, init = "last.par.best",
              seed = 101,
              lower = (1-0.4*sign(lpb))*lpb,
              upper = (1+0.4*sign(lpb))*lpb)
              
                  

## could also code this directly in Stan (https://mc-stan.org/docs/stan-users-guide/odes.html;
## also https://mpopov.com/tutorials/ode-stan-r/, https://shug3502.github.io/blog/DifferentialEqnsStan, ...)
