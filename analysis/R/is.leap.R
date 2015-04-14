is.leap <- function(year=c()) {

  # argument: year (YYYY)
  # returns: leap year T/F (logical)
  
  is.leap <- ((((year %% 4) == 0) & ((year %% 100) != 0)) | ((year %% 400) == 0))
  
  return(is.leap)
}
