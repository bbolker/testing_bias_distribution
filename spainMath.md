
The testing hazard for negatives is $b$, corresponding to a probability $p_n = 1-B$; for positives it is $b+h$, corresponding to a probability of $p_p = 1-BH$.

Let's say that group $k$ has size $n_k$ of whom $t_k$ are tested, and has hazard multiplier $\phi$ (so that $\phi_p = H$ and $\phi_n = 1$.

The probability of $t_k$ tests, given the parameters, is

$$ P_k = \binom{n_k}{t_k} (1-\phi_kB)^{t_k} (\phi_kB)^{n_k-t_k}$$

$$ \propto (1-\phi_kB)^{t_k} B^{n_k-t_k} \equiv L_k,$$

after dropping factors that don't depend on $B$.

The overall likelihood is:

$$L = \prod_k{(1-\phi_kB)^{t_k} B^{n_k-t_k}}$$

$$ =  \prod_k{(1-\phi_kB)^{t_k} B^{\sum(n_k-t_k)}}$$

$$ =  \prod_k{(1-\phi_kB)^{t_k} B^{N-T}}$$

Thus, the maximum-likelihood estimate of $B$ does not depend on the group sizes $n_k$, but only their sum the population size. This is because in this parameterization, an increase in $B$ makes each negative test more likely in the same proportion, thus has the same relative effect on their product no matter how they are distributed. Note that the $n_k$ \emph{do} affect the original probability, via the $\prod \phi_k^{n_k}$ that we scaled out.

We can then write log likelihood

$$ \ell = \sum_k t_k \log(1-\phi_kB) + (N-T) \log(B)$$

and look for zeros of 

$$ \ell_B = - \sum_k \phi_k \frac{t_k}{1-\phi_kB} + (N-T)/B$$

In terms of our original two groups, we could write this for example as:

$$ \ell_B = - \frac{x}{1-B} - H\frac{y}{1-HB} + (N-x-y)/B,$$

where $x$; $y$ is the number of negative and positive tests observed.

We have solved this equation in computer algebra. The numerator is quadratic, with the larger solution corresponding to a likelihood minimum (and is presumably $>1$ for sensible parameters).
