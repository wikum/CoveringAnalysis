

tryCatch({
  cat("\n PAIR NETWORK \n ===================== \n")
  source("1.R")
}, error = function(e){print(e)})

tryCatch({
  cat("\n SAMPLE NETWORKS \n ===================== \n")
  source("2.R")
}, error = function(e){print(e)})



