
We made some test data and tried to fit with nlme::gls(), but we don't believe it's actually using the formula (e.g., in weights=~x) as a formula, just as a fixed variance for weighting.

We are also struggling to understand why bs() is producing really weird results compared to ns, including slope discontinuities beyond what seems reasonable, we should learn more about this.
