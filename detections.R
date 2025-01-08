library(shellpipes)
manageConflicts()

library(purrr)
library(readr)
library(dplyr)

dl <- csvReadList(na = c("", "NA", "N.R"))
names(dl) <- sub(".*_[12][0-9]", "FY", names(dl))
names(dl) <- sub("/.*", "", names(dl))
print(names(dl))

warn <- (map(dl, problems)
	|> bind_rows(.id="file")
)
print(warn)

nl <- map(dl, names)

fields <- (nl |> unname() |> unlist() |> unique())
testFields <- grep("_tests", fields, value=TRUE)
confFields <- grep("_positive_tests", testFields, value=TRUE)
testFields <- setdiff(testFields, confFields)

tsvSave(arrange(data.frame(confFields)))
print(testFields)
 
