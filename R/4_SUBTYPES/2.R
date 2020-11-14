
# aberration probabilities

sink("2.log.txt", split=TRUE)

rm(list=ls())
gc()

tryCatch({
  
  set.seed(1)
  
  l1 = load("../OBJ/cover.source.rda")
  l2 = load("../DATA/TCGA/pheno.data.rda")
  l3 = load("../DATA/SOURCE/mat.list.source.rda")
  
  subtypes.list = list(
    list(tissue = "Breast",
         subtype = "PAM50",
         classes = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like")
    ),
    list(tissue = "Lung",
         subtype = "smoking",
         classes = c("Smoker", "RecentlyReformed", "Reformed", "NonSmoker")
    ),
    list(tissue = "Colon",
         subtype = "CRIS",
         classes = c("CRISA", "CRISB", "CRISC", "CRISD", "CRISE")
    ),
    list(tissue = "Prostate",
         subtype = "primary",
         classes = c("3", "4", "5")
    )
  )
    
  for(i in 1:length(subtypes.list)){
    
    tissue = subtypes.list[[i]]$tissue
    subtype = subtypes.list[[i]]$subtype
    classes = subtypes.list[[i]]$classes
    
    mat = mat.list.source[[tissue]]
    sol = list.cover.source[[tissue]]$core
    
    pheno = pheno.list[[tissue]]
    
    all(rownames(pheno) == colnames(mat))  
    
    groups = as.character(pheno[, subtype])
    
    mat.probs = sapply(classes, function(x){
      sel = which(groups == x)
      z = mat[sol, sel]
      round(rowMeans(z), 3)
    })
    
    print(mat.probs)
    
  }
  

}, error = function(e){ print(e) })

sink()








