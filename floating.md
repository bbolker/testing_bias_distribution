## 2025 Dec 13 (Sat)

Make some sort of simulator, and then move it to tmb world?

Thinking right now that macpan is not going to be perfect for the floating baseline approach, maybe start a parallel initiative with raw sims and bbmle

## Prompt

I am running a discrete-time simulator in R. My state is a named list of state variables. What is an efficient way to accumulate the list into a data frame?

I want to pass the list for clean inner code, but also use vector updating. For example, I now have a function that does

	grad <- c(1, rateFun(states, params))
	return(unlist(states) + grad*params$deltat)

How can I efficiently get this series of states into a data frame? I don't mind if I accumulate in another form and turn to data frame at the end.
