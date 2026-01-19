
The testing hazard for negatives is $b$, corresponding to a probability $p_n = 1-B$; for positives it is $b+h$, corresponding to a probability of $p_p = 1-BH$.

The probability that $T$ people are tested in a group is proportional to $p^T (1-p)^{(X-T)}$, where $X$ is the size of the group. The log-likelihood then satisfies
$$l \propto T \log(p)+(X-T)\log(1-p)$$
So with the current prevalence $Y$ for the number of observed positive tests $pos$,  the corresponding log-likelihood satisfies:$$l_p \propto pos\times \log(1-BH)+(Y-pos)\times \log(BH)$$
Similarly, for the number of observed negative tests $neg$, corresponding log-likelihood satisfies:$$l_n \propto neg\times \log(1-B)+(N-Y-neg)\times \log(B)$$
The baseline hazard $b$ that maximize the log-likelihood should satisfies the equation:$$\frac{d}{dB}(l_p+l_n)=0$$
so we have $$\begin{split}
(-\frac{pos}{1-BH}+\frac{Y-pos}{BH})H+(-\frac{neg}{1-B}+\frac{N-Y-neg}{B})&=0
\\
\Leftrightarrow  \frac{-B^2 N H + ((N - neg)H + N - pos)B - N + neg + pos}{B(1-B)(1-BH)}&=0
\\
\Leftrightarrow  N H B^2 - ((N - neg)H + N - pos)B +(N - neg - pos) &=0
\end{split}$$
As a result, we have two potential $B$ solution that maximize the likelihood $$B=\frac{1}{2NH}\Big[(N-neg)H+N-pos \pm \sqrt{((N-neg)H+N-pos)^2-4N H (N-pos-neg)} \Big]$$
Currently both solution seems to be positive (feasible) and be independent of $Y$!
