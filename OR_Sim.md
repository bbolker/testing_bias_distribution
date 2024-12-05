If we don't allow for false readings, then the number of positive tests is just the number of positive people tested (and the number of negative tests is just the number of negative people tested). So the most straightforward way to calculate the likelihood is dbinom(posTests, posPeople, posTestingProb)*dbinom(negTests, negPeople, negTestingProb).

posPeople and negPeople come from the current sim, the Probs come from the parameters, and the Tests come from the “data”. so we don't need to simulate tests in the current sim.

This clarifies also that the most likely way for the model to fall off the cliff would be if posTests (from the “data”) exceeds posPeople (from the current sim).
