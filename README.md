# testing_bias_distribution

Models to understand dynamics of test positivity under epidemic dynamics and biased testing

-   [web version](https://bbolker.github.io/testing_bias_distribution/)
-   [beta notes](docs/testing_distrib.html)
-   [simple interaction model](simple.md)
-   [RTMBode/fitode example](fitodeMCMC.R): not directly relevant, but hopefully a useful example of how to do trajectory-matching in `fitode` or `RTMBode`

Related projects. It's not super-clear we have any of sufficient interest, but JD is willing to share. We have a mathy publication with an earlier post-doc (Ali Gharouni, nice guy, now at PHO); Mike Li did some stuff during Covid which probably had some ideas but will be hard to figure out; JD advised some South African students who have a public repo, but it's not clear how far they got.

## files

-   `corrCheck.R`: numerical calculations of correlation showing non-monotonicity
-   `Expected_Test_positivity_figure.R`: test positivity as a function of testing focus/dispersion \phi, testing proportion t, prevalence i
-   `simple.md`: basic notes on multiplex testing and cross-effects of diseases on each others' positivity/number of positive tests
-   `testing_funs.R`: basic machinery for computing expected positivity based on the Beta model
-   `testing_distrib.rmd`: exploration of properties of `testing_funs.R`
-   `Logspace_comparing_methods.R`: the file used to generate the pictures for log-space comparing.
-   `Qbeta_Issues.R`: this is created regarding the Qbeta & qbeta issues.
-   `Single_variable_SIR.R`: a numerical simulation of simple version SIR trajectory.
-   `inc_testing_positivity-ratio.R`: As discussed, this file is used to generate figures for ratio of inc/test_positivity_proportion under different combination of inc, phi and testing_proportion.
-   `Average_i-phi_idea.R`: Unsuccessful experiments. I'll try something else but please ignore this for now.

### In Pix
- `ratio_inc-pos_prop.png`: heat map of whole parameter space for inc and phi, grouped by different testing proportion. Color is the log-ratio. For the lwr=1 area, the result is replaced by NaN since the they are clearly wrong and also for better illustration of the heatmap color.
- `ratio_focus_inc_pos_prop.png`: zoomed version to a more reasonable inc range (0,0.25)
- `ratio_phi_slice.png`: Slice of the heat map for different phi (color), grouped by test proportions. Y-axis is ratio, X-axis is inc. For better illustration of changing in phi.
ratio_inc_slice.png: Slice of the heat map for different inc (group), colored by test proportions. Y-axis is ratio, X-axis is phi. For better illustration of changing in inc.

-   `ratio_inc-pos_prop.png`: heat map of whole parameter space for inc and phi, grouped by different testing proportion. Color is the log-ratio. For the lwr=1 area, the result is replaced by NaN since the they are clearly wrong and also for better illustration of the heatmap color.
-   `ratio_focus_inc_pos_prop.png`: zoomed version to a more reasonable inc range (0,0.25)
-   `ratio_phi_slice.png`: Slice of the heatmap for different phi (color), grouped by test proportions. Y-axis is ratio, X-axis is inc. For better illustration of changing in phi. ratio_inc_slice.png: Slice of the heatmap for different inc (group), colored by test proportions. Y-axis is ratio, X-axis is phi. For better illustration of changing in inc.
