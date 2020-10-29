##### help functions 

#### compute divergent probability for list of matrix
DivProb <- function(Mat){
  temp <- lapply(Mat, function(x){
    x <- abs(x)
    temp1 <- apply(x, 1, function(y){
      sum(y==1)/length(y)
    })
    return(temp1)
  })
  names(temp) <- names(Mat)
  return(temp)
}


### build motif with module, where source module g  contains all pairs start with source gene g
### target module contains all pairs start with target gene, 
### Mat: list of binary expression for pairs across all tissues 
MotifExp <- function(Module, Mat){
  
  temp <- lapply(1:length(Module), function(i){
    temp <-Module[[i]]
    mat <- Mat[[i]]
    module_exp<-lapply(temp, function(x){
      if (length(x) <2){
        rr <- matrix(mat[x,], nrow=1)
      }else{
        rr <- mat[x,]
      }
      result <- apply(rr,2, max)
      return(result)
    })
    module_exp <- do.call(rbind, module_exp)
    colnames(module_exp) <- colnames(mat)
    return(module_exp)
    
  })
  
  names(temp) <- names(Mat)
  return(temp)
}


#####search for covering with alpha = 0, J = 1, 2,3, mat: list of matrix from pair, source, target level 
MinimalSigGene2 <- function(mat, alpha, J){
  
  ### alpha is the fraction of samples not covered by the signature
  Result_Motif<-lapply(1:length(mat), function(i){
    
    alpha <-alpha
    A <- mat[[i]]
    Asum <- sum(apply(A, 2, sum) <J)
    if (Asum <= alpha*ncol(A)){
      alpha=alpha
    }else{
      #### change since 2% experiment
      alpha <- Asum/ncol(A)
    }
    print(alpha)
    AA <- rbind(A, J*diag(ncol(A)))
    AA <- cbind(AA,c(rep(0,nrow(A)),rep(1, ncol(A))))
    
    library(gurobi)
    model <- list()
    model$A <- t(AA)
    model$obj <- c(rep(1, nrow(A)), rep(0,ncol(A)))
    model$modelsense <- 'min'
    model$modelname <- 'poolserch'
    model$rhs <- c(rep(J, ncol(A)),alpha*ncol(A))
    
    model$sense <- c(rep('>=', ncol(A)),'<=')
    model$vtype <- 'B'
    params <- list()
    params$PoolGap <- 0
    params$PoolSolutions <- 100000
    params$PoolSearchMode <-2
    
    result <- gurobi(model, params = params)
    
    ### select indicator for genes
    result_pool1 <- lapply(result$pool, function(x){x$xn[1:nrow(A)]})
    result_pool2 <- lapply(result$pool, function(x){x$xn[(nrow(A)+1):(nrow(A)+ncol(A))]})
    
    gene1 <- matrix(unlist(result_pool1), nrow = length(result_pool1), byrow = TRUE)
    sample1 <- matrix(unlist(result_pool2), nrow = length(result_pool2), byrow = TRUE)
    
    colnames(gene1) <- rownames(A)
    colnames(sample1) <- colnames(A)
    genesig_size <- result$objval
    genesigs_all <- nrow(gene1)
    return(list(genesig_size=genesig_size, gene1=gene1,sample1=sample1, alpha=alpha))
    
  })
  
  names(Result_Motif) <- names(mat)
  return(Result_Motif)
  
}






