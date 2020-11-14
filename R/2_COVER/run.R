

tryCatch({
  cat("\n PAIR COVERINGS \n ===================== \n")
  source("1.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n SOURCE COVERINGS \n ===================== \n")
  source("2.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n TARGET COVERINGS \n ===================== \n")
  source("3.R")
}, error = function(e){print(e)})


