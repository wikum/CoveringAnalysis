
# pair coverings

rm(list=ls())
gc()

require(lpcover)

l1 = load("../DATA/PAIR/mat.list.pair.rda")

mat.list.pair = mat.list.pair[-4]
names(mat.list.pair)[3] = "Kidney"

J = 1

list.cover.pair = lapply(mat.list.pair, function(mat){

  alpha = 1 - sum(colSums(mat) >= J)/ncol(mat)

  R = computeMinimalCovering(mat, alpha=alpha, maxsol=10, J=J, solver="gurobi")

  mat.sol = R$sol

  union = unique(as.vector(mat.sol))
  pvals = rowMeans(mat[union, ])
  sel.sol = which.max(apply(mat.sol, 2, function(x) sum(pvals[x]) ))[1]
  sol = mat.sol[, sel.sol]
  core = union[which(rowMeans(apply(mat.sol, 2, function(x) union %in% x )) == 1)]

  alpha = sum(colSums(mat[sol, ]) >= J)/ncol(mat)
  
  list(mat.sol=mat.sol,
       union=union,
       core=core,
       sol=sol,
       alpha=alpha)
  
})

save(list.cover.pair, file="../OBJ/cover.pair.rda")

print(
t(sapply(list.cover.pair, function(x) c(size=length(x$sol),
                                 union=length(x$union),
                                 core=length(x$core),
                                 coverage=round(x$alpha, 3)) 
       ))
)





