To prove the concept that we can gain something by fitting both series, one approach would be to try to fit an epidemic without fitting a baseline.

I had argued previously that fitting hazards is hard, but now that I'm actually working on practical approaches, I'm thinking that the floating baseline hazard may make things more stable; we're only fitting a single relative hazard. It's true that this is constrained: we need r+b>0 , where r is the relative hazard and b is the floating hazard. But in practice if we take b as the hazard for negative people, then we expect both b and r to stay comfortably positive.

2025 Dec 16 (Tue) I have sort of proved the concept for what I think is a silly example (using a risk ratio rather than an odds ratio). It looks like you can in fact get a reasonable amount of information by combining time series in a perfect world, at least. I'll make some notes about this (you can look at the sevilla/ subdirectory).

Now I'm trying to ask whether we can use math and uniroot to find the best hazard value for a set of observations under the assumption of fixed relative hazard. mathHazard.R … yes. So this means we can try to do a proof of concept for relative hazard.
