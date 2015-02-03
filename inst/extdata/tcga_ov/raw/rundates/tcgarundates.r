# TCGA Run dates script 
# January 28, 2015
# setwd("~/Documents/RTCGAToolbox/")
# Function to format date times
fdate <- function(dat){
  dat <- gsub("([0-9]{4})([0-9]{2})", "\\1-\\2-", dat)
  return(dat)
}

# all.dates <- rev(fdate(getFirehoseRunningDates()))
# an.dates <- rev(fdate(getFirehoseAnalyzeDates()))
# Saving last 5 dates and today's date (checked)
# dates <- data.frame(run = tail(all.dates, 5), analyze = tail(an.dates, 5), 
#    checked = as.character(rep(Sys.Date(),5)), stringsAsFactors = FALSE)

load("./rundates/rundates.rda")
if(dates$checked[nrow(dates)] != as.character(Sys.Date())) {
  dates <- rbind(dates, c(fdate(getFirehoseRunningDates(last = 1)), 
                 fdate(getFirehoseAnalyzeDates(last = 1)), as.character(Sys.Date())) )
  save(dates, file = "rundates.rda")
}

