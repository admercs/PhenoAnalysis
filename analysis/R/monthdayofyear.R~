dayofyear <- function(year=c(), month=c(), day=c()) {

  # arguments:
  # year : YYYY
  # month : MM
  # day : DD
  
  # returns:
  # day of year

  source("is.leap.R")
  
  if (is.leap(year)) {
    monthdays <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  } else {
    monthdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  }

  if (month > 1) {
    dayofyear <- sum(monthdays[1:month-1]) + day
  } else {
    dayofyear <- day
  }
  
  return(dayofyear)
}
