For now, we are considering a compartmental SIR like disease.

Currently, for both Beta distribution and Odd ratio model, at each time point $t$, the proportion of positive test $p$ and the proportion of all test $T$ are both determined by prevalence $Y=I/N$, which is the proportion of infected population at time point $t$.

$$ T= (1-Y) T_B + Y T_Y $$

$$ p= \frac{Y T_Y}{T}= \frac{1}{1+(\frac{1}{Y}-1)\frac{T_B}{T_Y}}$$

Based on $T$ and $p$, we determine the observed number of test $OT$ and observed number of positive test $OP$ at time $t$.

$$OT \sim \text{Binom}(N,T)$$

$$OP \sim \text{Binom}(OT,p)$$
Since $OT$ and $OP$ are both observed "incidence" at time $t$, a question would be if it is more reasonable to use incidence of infection $\frac{S(t)-S(t-1)}{N}$ when calculating the rates. 

What we are doing now with prevalence, is assuming that at each time point $t$, all active infected cases at the moment $I(t)$ are tested based on rate $T_Y$ and all other population $N-I(t)$ are tested based on baseline rate $T_B$. 

### Thoughts:
This might be more justifiable than using incidence, since being symptomatic/actively infected is a key factor for $T_Y \gg T_B$. From this perspective, it is problematic to directly link new infections to new positive tests (or with a delay?):
- Newly infected cases might be asymptomatic thus less possible to get tested
- While long infection duration time could also have more mild/ignorable symptomatic with less testing possibility due to self-recovery.). 
- Other reasons driven positive testing rate could be complicated and independent with time.
A harmonized uniform testing rate for all active infected population at the moment might be more reasonable.

With prevalence, the tendency of $OP$ is synchronized with prevalence $Y(t)$ and so does $OT$ in general as we assume $T_Y > T_B$.

$$ \frac{dT}{dt}=\frac{dT}{dY} \frac{dY}{dt}= (-T_B +T_Y)\times \frac{dY}{dt}$$

$$ \frac{dp}{dt}=\frac{dp}{dY}\frac{dY}{dt}= \frac{\frac{T_B}{T_Y}}{((Y-1)\frac{T_B}{T_Y}-Y)^2}\times\frac{dY}{dt}$$

However, intuitively, a large new infection should somehow boost the test number $OT$ and positive test $OP$ number (within some duration in the future). The peak of incidence might not overlap with peak of prevalence, due to recovery process. A draw back using prevalence is that we might not be able to catch this.

I think the key difference would be if the new infection dominantly derive the positive testing or the positive tests is more a mixture result of new and cumulative infection.
### Positive tests' impact on dynamic 
A related question would be if we need to consider impact to infection dynamic of those observed positive tests: 
- If they are positive in test, will they get recovered faster (take medication, treated, hospitalized, etc.) 
- And/or get less infective (self/administrative quarantined, wear masks, extra measurement etc.). 
- This impact might not be significant if $T_{prop}$ and/or $pos$ are small, thus can be ignored for simplicity.