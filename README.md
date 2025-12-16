# testing_bias_distribution

Models to understand dynamics of test positivity under epidemic dynamics and biased testing

-   [web version](https://bbolker.github.io/testing_bias_distribution/)
-   [beta notes](docs/testing_distrib.html)
-   [simple interaction model](simple.md)
-   [RTMBode/fitode example](fitodeMCMC.R): not directly relevant, but hopefully a useful example of how to do trajectory-matching in `fitode` or `RTMBode`

Related projects. It's not super-clear we have any of sufficient interest, but JD is willing to share. We have a mathy publication with an earlier post-doc (Ali Gharouni, nice guy, now at PHO); Mike Li did some stuff during Covid which probably had some ideas but will be hard to figure out; JD advised some South African students who have a public repo, but it's not clear how far they got.

## files
-   [`Notes_Nov10.Rmd`](Notes_Nov10.Rmd): @RichardSichengZhao 's summary notes for progress
-   [`corrCheck.R`](corrCheck.R): numerical calculations of correlation showing non-monotonicity
-   [`Expected_Test_positivity_figure.R`](Expected_Test_positivity_figure.R): beta model's test positivity as a function of testing focus/dispersion $\phi$, testing proportion t, prevalence i
-   [`simple.md`](simple.md): basic notes on multiplex testing and cross-effects of diseases on each others' positivity/number of positive tests
-   [`testing_funs.R`](testing_funs.R): basic machinery for computing expected positivity based on the Beta model and Hazard Ratio Model
-   [`testing_distrib.rmd`](testing_distrib.rmd): description and exploration of properties of `testing_funs.R`
-   [`Logspace_comparing_methods.R`](Logspace_comparing_methods.R): the file used to generate the pictures for log-space comparing 4 different methods of calculating positivity for Beta model.
-   [`Qbeta_Issues.R`](Qbeta_Issues.R): investigation of the Qbeta & qbeta issues.
-   [`Single_variable_SIR.R`](Single_variable_SIR.R): a numerical simulation of simple version SIR trajectory (fitODE later).
-   [`inc_testing_positivity-ratio.R`](inc_testing_positivity-ratio.R): As discussed, this file is used to generate figures for ratio of inc/test_positivity_proportion under different combination of inc, phi and testing_proportion.
-   [`Average_i-phi_idea.R`](Average_i-phi_idea.R): RZhao 's Unsuccessful experiments. Will try something else but please ignore this for now.
-   [`HR_Test_positivity_figure`](HR_Test_positivity_figure): Hazard Ratio model's test positivity curves as a function of testing Risk offsets \Phi, testing proportion T, prevalence V

* conceptNotes.md from JD in Sevilla

### In Pix
(TO DO)
