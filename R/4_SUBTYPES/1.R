
# motif subtype analysis 

rm(list=ls())
gc()

sink("1.log.txt", split=TRUE)

tryCatch({
  
  set.seed(1)
  
  library(rutils)
  library(lpcover)
  
  l1 = load("../DATA/PAIR/mat.list.pair.rda")
  mat.list.pair = mat.list.pair[-4]
  names(mat.list.pair)[3] = "Kidney"
  
  l2 = load("../OBJ/cover.pair.rda")
  
  l3 = load("../DATA/TCGA/pheno.data.rda")  
  pheno.list = pheno.list[-4]
  names(pheno.list)[3] = "Kidney"
  
  alphas = 0.1
  names(alphas) = paste("alpha", alphas, sep="")
  
  optim.subtypes.alpha.list = list()
  
  run_single_alpha = function(mat, groups, classes, alpha, maxsol, J){
    
    tb = table(groups)
    tb = tb[tb > 0]
    
    min.count = min(tb)
    
    sel.samples.list = lapply(names(tb), function(x) sample(which(groups == x), min.count))
    names(sel.samples.list) = names(tb)
    sel.samples = Reduce(c, sel.samples.list)
    
    out = lpcover::computeMinimalCovering(mat=mat[, sel.samples], alpha=alpha, maxsol=maxsol, J=J)
    out$sample.ids = colnames(mat[, sel.samples])
    out$alpha = alpha
    out$groups = groups[sel.samples]
    out$classes = classes
    
    sol = out$sol[, 1]
    out$prop = sapply(sel.samples.list, function(x) mean(colSums(mat[sol, x, drop=FALSE]) >= J) )
    
    out
    
  }
  
  N = 20
  
  # ----------------------------------------------------------------------------------
  # breast
  optim.subtypes.alpha.list$Breast = list()
  
  mat = mat.list.pair$Breast
  P = pheno.list$Breast
  all(colnames(mat) == rownames(P))
  
  # PAM50
  classes = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like")
  sel = which(P$PAM50 %in% classes & colSums(mat) > 0)
  groups = P$PAM50[sel]
  
  optim.subtypes.alpha.list$Breast$PAM50 = lapply(alphas, function(a){
    
    lapply(1:N, function(j){
      run_single_alpha(mat=mat[, sel], groups=groups, classes=classes, alpha=a, maxsol=1, J=1)
    })

  })

  rm(mat, P)
  
  # ----------------------------------------------------------------------------------
  
  # colon
  optim.subtypes.alpha.list$Colon = list()
  
  mat = mat.list.pair$Colon
  P = pheno.list$Colon
  all(colnames(mat) == rownames(P))
  
  # CRIS
  
  classes = c("CRISA", "CRISB", "CRISC", "CRISD", "CRISE")
  sel = which(P$CRIS %in% classes & colSums(mat) > 0)
  groups = P$CRIS[sel]
  
  optim.subtypes.alpha.list$Colon$CRIS = lapply(alphas, function(a){
    
    lapply(1:N, function(j){
      run_single_alpha(mat=mat[, sel], groups=groups, classes=classes, alpha=a, maxsol=1, J=1)
    })
    
  })
  
  rm(mat, P)
  
  # ----------------------------------------------------------------------------------
  
  # lung
  optim.subtypes.alpha.list$Lung = list()
  
  mat = mat.list.pair$Lung
  P = pheno.list$Lung
  all(colnames(mat) == rownames(P))

  # smoking
  
  classes = c("RecentlyReformed", "Reformed", "Smoker", "NonSmoker")
  sel = which(P$smoking %in% classes & colSums(mat) > 0)
  groups = P$smoking[sel]
  
  optim.subtypes.alpha.list$Lung$smoking = lapply(alphas, function(a){
    
    lapply(1:N, function(j){
      run_single_alpha(mat=mat[, sel], groups=groups, classes=classes, alpha=a, maxsol=1, J=1)
    })
    
  })
  
  
  rm(mat, P)
  
  # ----------------------------------------------------------------------------------
  
  # prostate
  optim.subtypes.alpha.list$Prostate = list()
  
  mat = mat.list.pair$Prostate
  P = pheno.list$Prostate
  all(colnames(mat) == rownames(P))
  
  # primary
  classes = 3:5
  sel = which(P$primary %in% classes & colSums(mat) > 0)
  groups = P$primary[sel]
  
  optim.subtypes.alpha.list$Prostate$primary = lapply(alphas, function(a){
    
    lapply(1:N, function(j){
      run_single_alpha(mat=mat[, sel], groups=groups, classes=classes, alpha=a, maxsol=1, J=1)
    })
    
  })
  
  rm(mat, P)
  
  # ----------------------------------------------------------------------------------
  
  save(optim.subtypes.alpha.list, file="../OBJ/optim.subtypes.alpha.list.rda")
  
  # ----------------------------------------------------------------------------------
  
  boxplot(t(sapply(optim.subtypes.alpha.list$Breast$PAM50$alpha0.1, function(z) z$prop)))
  
  boxplot(t(sapply(optim.subtypes.alpha.list$Colon$CRIS$alpha0.1, function(z) z$prop)))
  
  boxplot(t(sapply(optim.subtypes.alpha.list$Lung$smoking$alpha0.1, function(z) z$prop)))
  
  boxplot(t(sapply(optim.subtypes.alpha.list$Prostate$primary$alpha0.1, function(z) z$prop)))
  

}, error = function(e){ print(e) })

sink()








