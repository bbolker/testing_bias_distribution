## code extracted from testing_distrib.rmd; RTMB works better

```r
mle_out<- mle2(c~dbinom(prob=prop_pos_test(plogis(logit_i_0)*exp(r*tau),
                                           t, exp(log_phi),
                                 debug=FALSE, phiscale = "unconstrained"),
              size=t*N),
              start=list(logit_i_0=qlogis(true_pars["i_0"]),
                         r=true_pars["r"], log_phi=log(true_pars["phi"]) ),
              data= list(tau=dd$tau, c=dd$c, t=dd$t, N=true_pars["N"]),
              control=list(maxit=1000)
)
mle_est0 <- coef(mle_out)
```


```{r pp,cache=TRUE, warning = FALSE}
pp <- profile(mle_out)
```

```{r orig-profile-pic, warning = FALSE}
ggdd <- as.data.frame(pp)
ggplot(ggdd,aes(focal,abs(z))) + geom_point() + geom_line() +
  facet_wrap(~param,scale="free_x")+
  scale_y_continuous(limits=c(NA,5))
```

```{r estimates}
invlink <- function(x) c(plogis(x[1]),x[2],exp(x[3]))
tt <- (broom::tidy(mle_out)
  %>% mutate(lwr=estimate-2*std.error,
             upr=estimate+2*std.error)
  %>% select(term,estimate,lwr,upr)
  %>% mutate_if(is.numeric,invlink)
  %>% mutate(term=c("i_0","r","phi"))
  %>% mutate(true=true_pars[c("i_0","r","phi")])
)
```
