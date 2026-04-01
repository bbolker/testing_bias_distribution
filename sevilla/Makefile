## This is a temporary subdirectory of testing_bias_distribution
## So JD can go a bit crazy on his Sevilla trip

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -cl "vmt"

######################################################################
## Look here 2026 Jan 07 (Wed)

## Trying to rebuild here with a binomial pathway, and hazard ratios.
## At some point move the good stuff to a new directory (floating/?)
## Or the main directory?
sirBin.Rout: env.R sirFuns.rda binFuns.rda
	$(pipeR)

## Working on the concern scenario on an airplane 2025 Dec 19 (Fri)
## base.sim.Rout: sim.R base.params.R binFuns.R
impmakeR += sim
%.sim.Rout: sim.R %.params.rda sirBin.rda
	$(pipeR)

## Flatten the parameters (constant concern)
## This should be chained
flat.params.Rout: flat.R base.params.rda
	$(pipeR)

## base.obs.Rout: obs.R base.params.R
## flat.obs.Rout: obs.R
%.obs.Rout: obs.R %.sim.rda %.params.rda
	$(pipeR)

## Plots (not implemented, maybe not needed, see .obs)
## base.simPlots.Rout: simPlots.R base.params.R
%.simPlots.Rout: simPlots.R %.sim.rda
	$(pipeR)

## Set starting points for fits
impmakeR += start
%.start.Rout: start.R %.params.rda
	$(pipeR)

######################################################################

## Fit number positive
## flat.posFit.Rout: posFit.R binFuns.R
## base.posFit.Rout: posFit.R binFuns.R
%.posFit.Rout: posFit.R sirBin.rda %.sim.rds %.start.rda
	$(pipeR)

## Fit observed prevalence NOT IMPLEMENTED
## flat.prevFit.Rout: prevFit.R binFuns.R

## Fit test results
## flat.testFit.Rout: testFit.R binFuns.R
## base.testFit.Rout: testFit.R binFuns.R
%.testFit.Rout: testFit.R sirBin.rda %.sim.rds %.start.rda
	$(pipeR)

######################################################################

## Emergency stuff for talk
## Now split between poisson-style (deprecated) and binomial style

sirPois.Rout: env.R sirFuns.rda poisFuns.rda

Sources += $(wildcard *.R)
autopipeR = TRUE

## sirFuns.Rout: sirFuns.R

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

Makefile: | makestuff
makestuff:
	ln -s ../$@ . 

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/pdfpages.mk

-include makestuff/git.mk
-include makestuff/visual.mk
