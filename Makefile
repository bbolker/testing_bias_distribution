## This is https://github.com/bbolker/testing_bias_distribution


current: target
-include target.mk
Ignore = target.mk

-include makestuff/perl.def

######################################################################
## Modularized float-baseline fitting 2026 Mar 08
## Moved from /sevilla on 2026 Apr 01

Sources += $(wildcard *.R)
autopipeR = TRUE

## Parameter setup
## const.params: fixed baseline
## float.params: time-varying baseline

## Generate observed data from parameters
## const.data.Rout: float.data.sim.R const.params.rda
## float.data.Rout: float.data.sim.R float.params.rda
%.data.Rout: float.data.sim.R %.params.rda
	$(pipeRcall)

%.pois.data.Rout: pois.data.sim.R %.params.rda
	$(pipeRcall)

## Hack around a small chaining problem Mar 2026
float_update: float.params.Rout float.data.Rout

## Time series curves of simulated data
## const.plot.Rout: plot.dataview.R const.data.rda
## float.plot.Rout: plot.dataview.R float.data.rda
## %.plot.Rout: plot.dataview.R %.data.rda
##  	$(pipeRcall)

%.plot.Rout: pois.dataview.R %.data.rda
	$(pipeRcall)

## Fixed fitting with non-varying B_lik.
## FIXME: It is now using true values as starting values for fitting, should test the performance on vary starting point
## FIXME: Fit float baseline data with fixed mechanism generate NaNs. No priority as its performance is expected to be bad
const.fixed.fit.Rout: fixed.macpan.fit.R const.data.rda
	$(pipeRcall)

const.fixed.pois.fit.Rout: fixed.pois.fit.R const.pois.data.rda
	$(pipeRcall)
## Flex fitting with varying B_lik
## FIXME: It is now using true values as starting values for fitting, should test the performance on vary starting point
%.flex.fit.Rout: flex.macpan.fit.R %.data.rda
	$(pipeRcall)

## Compare best fitted curve with data
## Labeling things better: do some points shapes and line types: flipping the dot and lines
%.check.Rout: check.fit.R %.fit.rda
	$(pipeRcall)

#### note the fixed vs flex, float and const
## float.flex.check.Rout: Float baseline data, Flexible fitting mechanism
## const.flex.check.Rout: Const baseline data, Flexible fitting mechanism
## const.fixed.check.Rout: Constant baseline data, Fixed fitting mechanism

######################################################################

######################################################################

## meet 2025 Dec 11 (Thu)

Ignore += pix/
## Richard and Ben
sir_seasonal_test.Rout: sir_seasonal_test.R sirs_seasonal/tmb.R
	$(pipeR)

Sources += floating.md


## Try to simulate with macpan
simTest.Rout: simTest.R
	$(pipeR)

## mpFitting.html: mpFitting.md

######################################################################

subdirs += sevilla

sirFuns.Rout: sirFuns.R
	$(pipeR)

## sirSimp.Rout: sirSimp.R sirFuns.rda
sirSimp.Rout: sirSimp.R sirFuns.rda
	$(pipeR)

sirSimp.tplot.Rout: testplot.R  sirSimp.rds
	$(pipeR)

######################################################################

## log-diff_cdf-log_simp.png
## log-diff_simp-log_simp.png
## Test_positivity_vs_test_proportion_phi_inc.png

## OR_Sim.Rout: OR_Sim.R OR_Sim.md
OR_Sim_SIR_MacPan.Rout: OR_Sim_SIR_MacPan.R
	$(pipeR)

## OR_Sim.Rout: OR_Sim.R OR_Sim.md
OR_Sim.Rout: OR_Sim.R
	$(pipeR)

######################################################################

Ignore += docs *.html
Sources += $(wildcard *.rmd *.qmd *.md)

testing_distrib.html: testing_distrib.rmd

%.html: %.rmd
	echo "rmarkdown::render(\"$<\")" | R --slave

%.html: %.md
	$(pandocs)

%.html: %.qmd
	quarto render $<

docs/%.html: %.html
	cp $< $@

######################################################################

Sources += betaParams.md simple.md README.md

######################################################################

Sources += $(wildcard *.R)

## Not working yet 2024 Sep 19 (Thu)
corrCheck.Rout: corrCheck.R
	$(pipeR)

Expected_Test_positivity_figure.Rout: Expected_Test_positivity_figure.R testing_funs.R

testing_distrib.html: testing_distrib.rmd testing_funs.R

##### 2024 Oct 16 (Wed)

Logspace_comparing_methods.Rout: Logspace_comparing_methods.R
inc-testing_positivity-ratio.Rout: inc-testing_positivity-ratio.R

######################################################################

## Robust binomial?

## drbinom.md
drbinom.Rout: drbinom.R
	$(pipeR)

######################################################################

### Odds ratios

Sources += or.md $(wildcard *.mac)

## or.mac.out: or.mac or.md
## or_mac.html: or_mac.rmd

or.mac.tex: or.mac.out mactex.pl
	$(PUSH)

######################################################################

Sources += spainReport.md

## sudo npm install -g markdown-cli-renderer

## Still having md-math problems 2026 Jan 19 (Mon)
spainMath.html: spainMath.md
	pandoc $< -s --mathjax \
	--css=https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/5.2.0/github-markdown.min.css \
	-o $@


spainMath.gfmview: spainMath.md
	grip $< 

## spainMath.pdf: spainMath.md
spainMath.pdf: spainMath.tex
	pdflatex $<
Ignore += spainMath.tex
spainMath.tex: spainMath.md
	$(pandocs)

Sources += $(wildcard *.max)
spainMath.maxima: spainMath.max

## rvdss stuff moved to ariCanada

## Some implicit curves of positivity vs. test proportion
## First attempt
orCurves.Rout: orCurves.R

## Now modularized a bit
## Original grid (two prevalences, three shapes)
orGrid.Rout: orGrid.R

## Try to converge on a point
## orConv.Rout: orConv.R

## orConv.compPlots.Rout: compPlots.R orConv.R
## orGrid.compPlots.Rout: compPlots.R 
%.compPlots.Rout: compPlots.R %.rds
	$(pipeR)

######################################################################

is.Rout: is.R
	$(pipeR)

## Calculate best hazard for fixed hazard ratio?
mathHazard.Rout: mathHazard.R
	$(pipeR)

## Beta illustrations
betaIllus.Rout: betaIllus.R

######################################################################

Ignore += $(subdirs)
alldirs += $(subdirs)

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/02.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone --depth 1 $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/pipeR.mk
-include makestuff/max.mk
-include makestuff/texj.mk
-include makestuff/pdfpages.mk

-include makestuff/git.mk
-include makestuff/gitbranch.mk
-include makestuff/visual.mk
