# testing_bias_distribution

Models to understand dynamics of test positivity under epidemic dynamics and biased testing

- [web version](https://bbolker.github.io/testing_bias_distribution/)
- [beta notes](docs/testing_distrib.html)
- [simple interaction model](simple.md)

Related projects. It's not super-clear we have any of sufficient interest, but JD is willing to share. We have a mathy publication with an earlier post-doc (Ali Gharouni, nice guy, now at PHO); Mike Li did some stuff during Covid which probably had some ideas but will be hard to figure out; JD advised some South African students who have a public repo, but it's not clear how far they got.

## files

- `corrCheck.R`: numerical calculations of correlation showing non-monotonicity
- `Expected_Test_positivity_figure.R`: test positivity as a function of phi, t, i
- `simple.md`: basic notes on multiplex testing and cross-effects of diseases on each others' positivity/number of positive tests
- `testing_funs.R`: basic machinery for computing expected positivity based on the Beta model
- `testing_distrib.rmd`: exploration of properties of `testing_funs.R`
