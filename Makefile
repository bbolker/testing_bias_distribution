## This is https://github.com/bbolker/testing_bias_distribution

current: target
-include target.mk
Ignore = target.mk

-include makestuff/perl.def

######################################################################

## log-diff_cdf-log_simp.png
## log-diff_simp-log_simp.png
## Test_positivity_vs_test_proportion_phi_inc.png

## OR_Sim-Fit.Rout: OR_Sim-Fit.R

## OR_Sim.Rout: OR_Sim.R OR_Sim.md

######################################################################

Ignore += docs *.html
Sources += $(wildcard *.rmd *.qmd *.md)

testing_distrib.html: testing_distrib.rmd

%.html: %.rmd
	echo "rmarkdown::render(\"$<\")" | R --slave

%.html: %.qmd
	quarto render $<

docs/%.html: %.html
	cp $< $@

######################################################################

Sources += betaParams.md simple.md README.md

######################################################################

## If people don't object to having make install makestuff, you can activate the 00.stamp line below and this stuff will work.
## Or we could make more simple rules like those above

autopipeR = defined

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

## Canada data repo; developing now for MD talk. Need to clean up at some point

Ignore += rvdss_canada
rvdss_canada.update: | rvdss_canada
	cd $| && git pull
rvdss_canada:
	git clone https://github.com/dajmcdon/rvdss-canada $@
rvdss_canada/data: | rvdss_canada ;

Ignore += detections.csv
## Can't easily be processed together bc pdm09 stuff changes (maybe other stuff as well)
detections.csv: rvdss_canada/data
	cat $</data/season*/respiratory_detections.csv > $@

## Just summarizing right now
## We have 12 years of two major things plus hcov; only 1.5 years of sars afaict
detections.Rout: detections.R $(wildcard rvdss_canada/data/season*/respiratory_detections.csv)
	$(pipeR)

## Playing around only (2024 flu year is over) 
lastYear.Rout: lastYear.R rvdss_canada/data/*_*_2024/respiratory_detections.csv
	$(pipeR)

lastYear.plots.Rout: firstplot.R lastYear.rds
	$(pipeR)

## For years with sarscov2 (flu years 2023 and 2024; 2025 has a different name)
## 2023.fourpath.Rout: fourpath.R
## Getting very hacky! Some years don't have the flu total...
impmakeR += fourpath
%.fourpath.Rout: fourpath.R rvdss_canada/data/*_*_%/respiratory_detections.csv
	$(pipeR)

impmakeR += twopath
## 2014.twopath.Rout: twopath.R
%.twopath.Rout: twopath.R rvdss_canada/data/*_*_%/respiratory_detections.csv
	$(pipeR)

impmakeR += firstplot
## 2024.fourpath.firstplot.Rout: firstplot.R
## 2023.fourpath.firstplot.Rout: firstplot.R
## 2014.twopath.firstplot.Rout: firstplot.R twopath.R
%.firstplot.Rout: firstplot.R %.rda
	$(pipeR)

impmakeR += secondplot
## 2024.fourpath.secondplot.Rout: secondplot.R
## 2023.fourpath.secondplot.Rout: secondplot.R
## 2014.twopath.secondplot.Rout: secondplot.R twopath.R
%.secondplot.Rout: secondplot.R %.rda
	$(pipeR)

rvdssYears = 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023 2024
rvdssBrin = $(rvdssYears:%=%.twopath.firstplot.Rout.pdf)
rvdssBrin.pdf: $(rvdssBrin)
	$(pdfcat)

brin.Rout: brin.R 2014.twopath.R

######################################################################

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
