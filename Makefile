## This is https://github.com/bbolker/testing_bias_distribution

current: target
-include target.mk
Ignore = target.mk

-include makestuff/perl.def

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

## If people don't object to having make install makestuff, you can activate the 00.stamp line below and this stuff will work.
## Or we could make more simple rules like those above

## autopipeR = defined

Sources += $(wildcard *.R)

## Not working yet 2024 Sep 19 (Thu)
corrCheck.Rout: corrCheck.R

Expected_Test_positivity_figure.Rout: Expected_Test_positivity_figure.R testing_funs.R

testing_distrib.html: testing_distrib.rmd testing_funs.R

##### 2024 Oct 16 (Wed)

Logspace_comparing_methods.Rout: Logspace_comparing_methods.R
inc-testing_positivity-ratio.Rout: inc-testing_positivity-ratio.R

######################################################################

### Odds ratios

Sources += or.md $(wildcard *.mac)

## or.mac.out: or.mac or.md
## or_mac.html: or_mac.rmd

or.mac.tex: or.mac.out mactex.pl
	$(PUSH)

######################################################################

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

## Beta illustrations

betaIllus.Rout: betaIllus.R

######################################################################

Ignore += $(subdirs)

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
-include makestuff/pdfpages.mk

-include makestuff/git.mk
-include makestuff/gitbranch.mk
-include makestuff/visual.mk
