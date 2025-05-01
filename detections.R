library(shellpipes)
manageConflicts()

library(purrr)
library(readr)
library(dplyr)

## Read in a list of files and make sure there are no readr “problems”
navalues <- c("", "NA", "N.R", "not tested", "not available")
dl <- csvReadList(na = navalues, show_col_types=FALSE)
warn <- (map(dl, problems)
	|> bind_rows(.id="file")
)
stopifnot(nrow(warn)==0)

## Give data _sets_ short names by flu year.
## A flu year goes from summer to summer is named after the later year
names(dl) <- sub(".*_[12][0-9]", "FY", names(dl))
names(dl) <- sub("/.*", "", names(dl))
print(names(dl))

## Collect and simplify names of data _fields_ using a weird corrections file
## Should we be using more standard “data dictionary” methods?
## Note that we are sequentially replacing, and that it seems a bit weird.
orig <- fields <- map(dl, names) |> unname() |> unlist() |> unique()
print(fields)
corr <- tsvRead(show_col_types=FALSE)
print(problems(corr))
for (r in 1:nrow(corr)){
	fields <- sub(corr[[r, "pat"]], corr[[r, "rep"]], fields)
}

tsvSave(data.frame(orig, fields) 
	|> filter(orig != fields)
	|> arrange(fields)
)

testFields <- grep("_tests", fields, value=TRUE)
## testFields <- sub("_tests", "", testFields)
posFields <- grep("_positive_tests", testFields, value=TRUE)
totFields <- setdiff(testFields, posFields)

posFields <- sub("_positive_tests", "", posFields)
pos <- posFields |> unique() |> sort()
tsvSave(data.frame(pos), ext="pos.Rout.tsv")

totFields <- sub("_tests", "", totFields)
tot <- totFields |> unique() |> sort()
tsvSave(data.frame(tot), ext="tot.Rout.tsv")

