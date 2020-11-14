
# target coverings

require(lpcover)

rm(list=ls())
gc()

l1 = load("../DATA/TARGET/mat.list.target.rda")

mat.list.target = mat.list.target[-4]
names(mat.list.target)[3] = "Kidney"

J = 3

list.cover.target = lapply(mat.list.target, function(mat){

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

save(list.cover.target, file="../OBJ/cover.target.rda")

print(
  t(sapply(list.cover.target, function(x) c(size=length(x$sol),
                                            union=length(x$union),
                                            core=length(x$core),
                                            coverage=round(x$alpha, 3)) 
  ))
)

core.union = unique(unlist(sapply(list.cover.target, function(x) x$core)))

mem.core = 1 * sapply(list.cover.target, function(x) core.union %in% x$core)
rownames(mem.core) = core.union

print(mem.core)


