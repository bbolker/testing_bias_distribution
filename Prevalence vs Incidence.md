Currently, for both Beta distribution and Odd ratio model, at each time point $t$, the proportion of positive test $pos$ and the proportion of all test $T_{prop}$ are both determined by prevalence $pY=I/N$, which is the proportion of infected population at time point $t$.

Based on $T_{prop}$ and $pos$, we determine the observed number of test $OT$ and observed number of positive test $OP$ at time $t$.
$$OT \sim \text{Binom}(N,T_{prop})$$
$$OP \sim \text{Binom}(OT,pos)$$
Since $OT$ and $OP$ are both observed "incidence" at time $t$, a question would be if it is more reasonable to use incidence of infection when calculating the rates. 

What we are doing now with prevalence, is assuming that at each time point $t$, all active infected cases at the moment $I(t)$ are tested based on rate $T_Y$ and all other population $N-I(t)$ are tested based on baseline rate $T_B$. This is justifiable since being symptomatic/actively infected could be a 


A related question would be if we need to consider impact to infection dynamic of those observed positive tests: if they are positive in test, will they get recovered faster (take medication, treated, hospitalized, etc.) and/or get less infective (self/administrative quarantined, wear masks, extra measurement etc.). But this impact might not be significant if $T_{prop}$ and/or $pos$ are small, thus can be ignored for simplicity.