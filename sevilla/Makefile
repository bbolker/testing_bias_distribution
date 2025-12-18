## This is a temporary subdirectory of testing_bias_distribution
## So JD can go a bit crazy on his Sevilla trip

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -cl "vmt"

######################################################################
## Trying to rebuild here with a binomial pathway, and hazard ratios.
## At some point move the good stuff to a new directory (floating/?)

sirBin.Rout: env.R sirFuns.rda binFuns.rda

## base.sim.Rout: sim.R base.params.R
%.sim.Rout: sim.R %.params.rda sirBin.rda
	$(pipeR)

######################################################################
## Emergency stuff for talk
## All broken because I split sirFuns to make poisFuns.R
## Ah, let's fix it.

sirPois.Rout: env.R sirFuns.rda poisFuns.rda

Sources += $(wildcard *.R)
autopipeR = TRUE

## test.Rout: test.R sirFuns.R
test.Rout: test.R sirPois.rda

## Very simple SIR with Poisson observation
impmakeR += sirdata
## base.sirdata.Rout: sirdata.R base.params.R
%.sirdata.Rout: sirdata.R %.params.rda sirPois.rda
	$(pipeR)

## Plot it
## base.obsplot.Rout: obsplot.R
%.obsplot.Rout: obsplot.R %.sirdata.rds
	$(pipeR)

## Fit it?
## base.sirfit.Rout: sirfit.R %.sirdata.rds sirPois.R
%.sirfit.Rout: sirfit.R %.sirdata.rds sirPois.rda %.params.rda
	$(pipeR)

## K, boom. That worked great.

######################################################################

## Add variation in testing propensity 
## varfit has no odds for I and then just divides the probability!
## Best case, and still not working as well as before.
impmakeR += vardata
## base.vardata.Rout: vardata.R base.params.R
%.vardata.Rout: vardata.R %.params.rda sirPois.rda
	$(pipeR)

## Plot the scenario
## bigpop.varplot.Rout: varplot.R bigpop.params.R
%.varplot.Rout: varplot.R %.vardata.rda
	$(pipeR)
## A nice first cut, lots of complexity

## Can we fit it with a floating baseline? Yes! 
## But this is a fairly cherry-picked case
## bigpop.varfit.Rout: varfit.R vardata.R sirPois.R
## bigpop.varfit.Rout: varfit.R vardata.R sirPois.R
%.varfit.Rout: varfit.R %.vardata.rds sirPois.rda %.params.rda
	$(pipeR)
impmakeR += varfit

## Do a mixed fit.
## bigpop.mixed.Rout: sirfit.R vardata.R sirPois.R
%.mixed.Rout: sirfit.R %.vardata.rds sirPois.rda %.params.rda
	$(pipeR)
impmakeR += mixed

## Plot a fitted model
## bigpop.varcomp.Rout: compPlot.R bigpop.params.R
## base.varcomp.Rout: compPlot.R 
%.varcomp.Rout: compPlot.R %.params.rda %.varfit.rds sirPois.rda
	$(pipeR)

## We can use the same code to plot the naive fit
## bigpop.mixcomp.Rout: compPlot.R
## base.mixcomp.Rout: compPlot.R
%.mixcomp.Rout: compPlot.R %.params.rda %.mixed.rds sirPois.rda
	$(pipeR)

## Or also try to make a naive prevalence fit?

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## ln -s ../makestuff . ## Do this first if you want a linked makestuff ##

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/pdfpages.mk

-include makestuff/git.mk
-include makestuff/visual.mk
