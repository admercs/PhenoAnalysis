dayofyear <- function(year=c(), month=c(), day=c()) {

  # arguments:
  # year : YYYY
  # month : MM
  
  # returns:
  # ndayofmonth : number of days in this month

  source("is.leap.R")
  
  if (is.leap(year)) {
    monthdays <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  } else {
    monthdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  }

  ndayofmonth <- monthdays[month]

  return(ndayofmonth)
}
