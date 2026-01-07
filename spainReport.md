
I wanted to summarize some of the work I did for my trip to Spain. There is [a subdirectory](sevilla/), and I kind of feel like it's a dead end. But there is some valuable work that we should move back to here.

My goal, which I'm hoping to make some progress on in the short term, is proof of concept for a floating baseline approach. The idea is to imagine that we don't really know what's driving the negatives to get tested, but that the positives are experiencing the same hazard, along with an additional hazard. Although I may remain a big fan of log odds, I wasn't able to see any obstacles to using a more mechanistic, hazard-based approach here, so that's my current suggestion.

The goal is to make a model a simulation model with something interesting happening to the baseline, and then:
* fit an SIR to true incidence
* fail to fit well to observed incidence and/or positivity
* fit well again by using both data streams.

There is some code that simulates a level of concern driven by both incidence and fatigue. I don't feel strongly about it, but I'm pointing at it in case it might be useful. You should look at the first substantive section of [the subdirectory Makefile](sevilla/Makefile); it might be useful to look at the default sim](sevilla/

I also did some math about fitting binomial distributions in the hazard framework, which I suspect will be useful. 

I struggled to fit, and I'm hoping that it will be possible to shoehorn this stuff into the fitting framework that you guys have already made some painful progress on.
