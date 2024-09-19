
The classic way to parameterize the beta is with “size” or “dispersion” parameters. We're not going to say “dispersion” because it's misleading: more size means more information means less variance. The standard one is k=a+b (don't know if that is the standard name). The Dushoff one is s=ab/(a+b).

This project has been using “variation” parameters that increase with variability. The classic one would correspond to γ = 1/(a+b) and the alternative one would be θ = 1/a + 1/b. The current version of the code is also using φ = 1 - exp(γ). This is a bit confusing, but seems monotonic.

We calculate the shape parameters from the variation parameters as a=p/γ, or a=1/θ(1-p).
