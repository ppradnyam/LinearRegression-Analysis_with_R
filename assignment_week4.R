plotfunc <- function ()
{
path <- paste(getwd(),"/","outcome-of-care-measures",".csv", sep="")
outcome <- read.csv(path, colClasses="character")
head(outcome)
outcome[,11] <- as.numeric(outcome[,11])
hist(outcome[,11])
}


plotfunc()
