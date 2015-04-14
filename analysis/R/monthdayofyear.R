monthdayofyear <- function(year=c(), doy=c()) {

  # arguments:
  # year : YYYY
  # doy : DDD
  
  # returns:
  # month of year and day of month

  source("is.leap.R")
  
  if (is.leap(year)) {
    monthdays <- c(31,29,31,30,31,30,31,31,30,31,30,31)
  } else {
    monthdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  }

  if (doy > sum(monthdays)) {
    stop("Day of Year is too large")
  }
  
  month <- 1
  while (sum(monthdays[1:month]) < doy) {
    month <- month + 1
  }

  if (month > 1 ) {
    day <- doy - sum(monthdays[1:(month-1)])
  } else {
    day <- doy
  }
  
  return(c(day,month))
}
