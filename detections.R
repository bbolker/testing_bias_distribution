library(shellpipes)
manageConflicts()

library(purrr)

dl <- csvReadList()
map(dl, summary)
