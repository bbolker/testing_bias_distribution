library(shellpipes)
manageConflicts()

library(purrr)
library(readr)
library(dplyr)

dl <- csvReadList(na = c("", "NA", "N.R"))
warn <- (map(dl, problems)
	|> bind_rows(.id="file")
)
print(warn)

## Names of data _sets_
names(dl) <- sub(".*_[12][0-9]", "FY", names(dl))
names(dl) <- sub("/.*", "", names(dl))
print(names(dl))

## Names of data _fields
fields <- map(dl, names) |> unname() |> unlist() |> unique()
corr <- tsvRead()
for (r in 1:nrow(corr)){
	fields <- sub(corr[[r, "pat"]], corr[[r, "rep"]], fields)
}

testFields <- grep("_tests", fields, value=TRUE)
confFields <- grep("_positive_tests", testFields, value=TRUE)
testFields <- setdiff(testFields, confFields)

tsvSave(arrange(data.frame(confFields)))
print(testFields)
 
