## This is https://github.com/bbolker/testing_bias_distribution

current: target
-include target.mk
Ignore = target.mk

######################################################################

## log-diff_cdf-log_simp.png
## log-diff_simp-log_simp.png
## Test_positivity_vs_test_proportion_phi_inc.png

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

Sources += betaParams.md simple.md README.md notes.md

######################################################################

## If people don't object to having make install makestuff, you can activate the 00.stamp line below and this stuff will work.
## Or we could make more simple rules like those above

autopipeR = defined

## Sources += $(wildcard *.R)

## Not working yet 2024 Sep 19 (Thu)
corrCheck.Rout: corrCheck.R

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

Makefile: makestuff/00.stamp
makestuff/%.stamp:
	- $(RM) makestuff/*.stamp
	(cd makestuff && $(MAKE) pull) || git clone --depth 1 $(msrepo)/makestuff
	touch $@

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
