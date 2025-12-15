## This is a temporary subdirectory of testing_bias_distribution
## So JD can go a bit crazy on his Sevilla trip

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -cl "vmt"

######################################################################

autopipeR = TRUE

## test.Rout: test.R sirFuns.R
test.Rout: test.R sirFuns.rda

## Very simple SIR with Poisson observation
impmakeR += sirdata
## base.sirdata.Rout: sirdata.R base.params.R
%.sirdata.Rout: sirdata.R %.params.rda sirFuns.rda
	$(pipeR)

## Plot it
## base.obsplot.Rout: obsplot.R
%.obsplot.Rout: obsplot.R %.sirdata.rds
	$(pipeR)

## Fit it?
## base.sirfit.Rout: sirfit.R %.sirdata.rds sirFuns.R
%.sirfit.Rout: sirfit.R %.sirdata.rds sirFuns.rda %.params.rda
	$(pipeR)

## K, boom. That worked great.

######################################################################

## Add variation in testing propensity 
impmakeR += vardata
## base.vardata.Rout: vardata.R base.params.R
%.vardata.Rout: vardata.R %.params.rda sirFuns.rda
	$(pipeR)

## Plot the scenario
## bigpop.varplot.Rout: varplot.R bigpop.params.R
%.varplot.Rout: varplot.R %.vardata.rds
	$(pipeR)
## A nice first cut, lots of complexity

## Can we fit it with a floating baseline? Yes! 
## But this is a fairly cherry-picked case
## bigpop.varfit.Rout: varfit.R vardata.R sirFuns.R
## bigpop.varfit.Rout: varfit.R vardata.R sirFuns.R
%.varfit.Rout: varfit.R %.vardata.rds sirFuns.rda %.params.rda
	$(pipeR)
impmakeR += varfit

## Do a mixed fit.
## bigpop.mixed.Rout: sirfit.R vardata.R sirFuns.R
%.mixed.Rout: sirfit.R %.vardata.rds sirFuns.rda %.params.rda
	$(pipeR)
impmakeR += mixed

## Plot a fitted model
## bigpop.varcomp.Rout: compPlot.R bigpop.params.R
## base.varcomp.Rout: compPlot.R 
%.varcomp.Rout: compPlot.R %.params.rda %.varfit.rds sirFuns.rda
	$(pipeR)

## We can use the same code to plot the naive fit
## bigpop.mixcomp.Rout: compPlot.R
## base.mixcomp.Rout: compPlot.R
%.mixcomp.Rout: compPlot.R %.params.rda %.mixed.rds sirFuns.rda
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
