

tryCatch({
  cat("\n SUBTYPE COVERAGE \n ===================== \n")
  source("1.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n SUBTYPE ABERRATION PROBABILITIES \n ===================== \n")
  source("2.R")
}, error = function(e){print(e)})



