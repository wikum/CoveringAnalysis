
#### use all normal to compute divergence
#### build binary expression for rna, dna, mut,cnv 

rm(list=ls())

### load libraries

library(SummarizedExperiment)
library(divergence)
library(knitr)
library(plyr)

source("util.R")

DATA_DIR = "../DATA/TCGA"

####load dataset computed from previous study 
load("~/Dropbox (MechPred)/GeneProject/ Gene project/TF/analysis/PathwayCommons/result3/TCGATernaryAll.rda")

###########################construct matrix for RNA, Mutation, CNV##############
#### read table 

BreastPheno <- read.table(sprintf("%s/BREAST/BRCA_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,dec = ".", header = TRUE, row.names = 1, check.names = FALSE)
Breast <- read.table(sprintf("%s/BREAST/HiSeqV2.gz", DATA_DIR), 
                     stringsAsFactors = FALSE, header = TRUE, row.names = 1)

ColonPheno <- read.table(sprintf("%s/COLON/COAD_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names = 1, header = TRUE)
Colon <- read.table(sprintf("%s/COLON/HiSeqV2.gz", DATA_DIR), 
                    stringsAsFactors = FALSE, header = TRUE, row.names = 1)

KidneyPheno <- read.table(sprintf("%s/KIDNEY/KIRC_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Kidney<- read.table(sprintf("%s/KIDNEY/HiSeqV2.gz", DATA_DIR), 
                    stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LiverPheno <- read.table(sprintf("%s/LIVER/LIHC_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Liver<- read.table(sprintf("%s/LIVER/HiSeqV2.gz", DATA_DIR), 
                   stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LungPheno <- read.table(sprintf("%s/LUNG/LUAD_clinicalMatrix.gz", DATA_DIR), 
                        sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Lung<- read.table(sprintf("%s/LUNG/HiSeqV2.gz", DATA_DIR), 
                  stringsAsFactors = FALSE, header = TRUE, row.names = 1)

ProstatePheno <- read.table(sprintf("%s/PROSTATE/PRAD_clinicalMatrix.gz", DATA_DIR), 
                            sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Prostate<- read.table(sprintf("%s/PROSTATE/HiSeqV2.gz", DATA_DIR), 
                      stringsAsFactors = FALSE, header = TRUE, row.names = 1)

##### combined different tissues and get interesection of samples between phenotype and RNAseq data

tissues <- c('Breast','Colon', 'Kidney',  'Liver', 'Lung', 'Prostate')
MultiRNAseq <- list(Breast, Colon, Kidney, Liver, Lung, Prostate)
MultiPheno <- list(BreastPheno, ColonPheno, KidneyPheno, LiverPheno, LungPheno, ProstatePheno)

#### change "-" to "."

MultiPheno <- lapply(MultiPheno, function(x){
  rownames(x) <- gsub('-','.',rownames(x))
  return(x)
})



for (i in 1:length(tissues)){
  
  temp <- intersect(colnames(MultiRNAseq[[i]]), rownames(MultiPheno[[i]]))
  MultiRNAseq[[i]] <- MultiRNAseq[[i]][, temp]
  MultiPheno[[i]]   <- MultiPheno[[i]][temp,]
  
}

names(MultiRNAseq) <- tissues
names(MultiPheno) <-  tissues


#### check RNAseq and Pheno are consistent 
lapply(1:length(tissues), function(i){
  all(rownames(MultiPheno[[i]])== colnames(MultiRNAseq[[i]]))
})




### compute divergence



MultiDiv <- lapply(1:length(tissues), function(i){
  
  Tempmat <- MultiRNAseq[[i]]
  TempPheno <- MultiPheno[[i]]
  
  norm_samp <- rownames(TempPheno)[TempPheno$sample_type=='Solid Tissue Normal']
  tumor_sample <- rownames(TempPheno)[TempPheno$sample_type=='Primary Tumor']
  
  baseMat <- as.matrix(Tempmat[, norm_samp])
  dataMat <- as.matrix(Tempmat[, tumor_sample])
  seMat.base = SummarizedExperiment(assays=list(data=baseMat))
  seMat = SummarizedExperiment(assays=list(data=dataMat))
  
  div <- computeUnivariateDigitization(seMat = seMat,
                                       seMat.base = seMat.base,
                                       parallel = TRUE,
                                       computeQuantiles = TRUE, 
                                       gamma = seq(1,99,5)/100,
                                       beta=0.95)
  
  
  return(div$Mat.div)
  
})


names(MultiDiv) <- tissues





#### Binary Expression for Mutation
#### read table 

BreastPheno <- read.table(sprintf("%s/BREAST/BRCA_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,dec = ".", header = TRUE, row.names = 1, check.names = FALSE)
Breast <- read.csv(sprintf("%s/BREAST/mutation_curated_wustl_gene.gz", DATA_DIR),
                   sep='\t',stringsAsFactors = FALSE, row.names = 1, header = TRUE)

ColonPheno <- read.table(sprintf("%s/COLON/COAD_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names = 1, header = TRUE)
Colon <- read.csv(sprintf("%s/COLON/mutation_bcm_gene.gz", DATA_DIR),
                  sep='\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

KidneyPheno <- read.table(sprintf("%s/KIDNEY/KIRC_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Kidney<- read.csv(sprintf("%s/KIDNEY/mutation_bcm_gene.gz", DATA_DIR),
                  sep='\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LiverPheno <- read.table(sprintf("%s/LIVER/LIHC_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Liver<- read.table(sprintf("%s/LIVER/mutation_bcm_gene.gz", DATA_DIR), 
                   sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LungPheno <- read.table(sprintf("%s/LUNG/LUAD_clinicalMatrix.gz", DATA_DIR), 
                        sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Lung<- read.table(sprintf("%s/LUNG/mutation_broad_gene.gz", DATA_DIR), 
                  sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)

ProstatePheno <- read.table(sprintf("%s/PROSTATE/PRAD_clinicalMatrix.gz", DATA_DIR), 
                            sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Prostate<- read.table(sprintf("%s/PROSTATE/mutation_broad_gene.gz", DATA_DIR), 
                      sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)


##### combined different tissues and get interesection of samples between phenotype and RNAseq data

tissues <- c('Breast','Colon', 'Kidney',  'Liver', 'Lung', 'Prostate')
MultiMat <- list(Breast, Colon, Kidney,  Liver, Lung, Prostate)
MultiPheno <- list(BreastPheno, ColonPheno, KidneyPheno,LiverPheno, LungPheno, ProstatePheno)

MultiPheno <- lapply(MultiPheno, function(x){
  rownames(x) <- gsub('-','.',rownames(x))
  return(x)
})

for (i in 1:length(tissues)){
  temp <- intersect(colnames(MultiMat[[i]]), rownames(MultiPheno[[i]]))
  MultiMat[[i]] <- MultiMat[[i]][, temp]
  MultiPheno[[i]]   <- MultiPheno[[i]][temp,]
  
}

names(MultiMat) <- tissues
names(MultiPheno) <-  tissues

lapply(1:length(tissues), function(i){
  all(rownames(MultiPheno[[i]])== colnames(MultiMat[[i]]))
})

#### check sample phenotype

lapply(MultiPheno, function(x){table(x$sample_type)})

### check unique values
lapply(MultiMat, function(x){unique(unlist(x))})

### get common samples with IDList
MultiMut <- lapply(1:length(MultiMat), function(i){
  Tempmat <- data.matrix(MultiMat[[i]])
  mut <- Tempmat
  ### remove samples with !0 and !1
  sampleid <- unique(which(mut !=0 & mut!=1, TRUE)[,2])
  
  if (length(sampleid) >0){
    
    return(mut[,-sampleid])
    
  }else{
    
    return(mut)
    
  }
})

names(MultiMut) <- tissues



#### ternary expression of the cnv 

#### read table 

BreastPheno <- read.table(sprintf("%s/BREAST/BRCA_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,dec = ".", header = TRUE, row.names = 1, check.names = FALSE)
Breast <- read.csv(sprintf("%s/BREAST/TCGA/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR),
                   sep='\t',stringsAsFactors = FALSE, row.names = 1, header = TRUE)

ColonPheno <- read.table(sprintf("%s/COLON/COAD_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names = 1, header = TRUE)
Colon <- read.csv(sprintf("%s/COLON/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR),
                  sep='\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

KidneyPheno <- read.table(sprintf("%s/KIDNEY/KIRC_clinicalMatrix.gz", DATA_DIR), 
                          sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Kidney<- read.csv(sprintf("%s/KIDNEY/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR),
                  sep='\t', stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LiverPheno <- read.table(sprintf("%s/LIVER/LIHC_clinicalMatrix.gz", DATA_DIR), 
                         sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Liver<- read.table(sprintf("%s/LIVER/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR), 
                   sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)

LungPheno <- read.table(sprintf("%s/LUNG/LUAD_clinicalMatrix.gz", DATA_DIR), 
                        sep='\t',stringsAsFactors = FALSE,row.names = 1, header = TRUE)
Lung<- read.table(sprintf("%s/LUNG/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR), 
                  sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)

ProstatePheno <- read.table(sprintf("%s/PROSTATE/PRAD_clinicalMatrix.gz", DATA_DIR), 
                            sep='\t',stringsAsFactors = FALSE, row.names=1,header = TRUE)
Prostate<- read.table(sprintf("%s/PROSTATE/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", DATA_DIR), 
                      sep='\t',stringsAsFactors = FALSE, header = TRUE, row.names = 1)


##### combined different tissues and get interesection of samples between phenotype and RNAseq data
tissues <- c('Breast','Colon', 'Kidney', 'Liver', 'Lung', 'Prostate')
MultiMat <- list(Breast, Colon,Kidney,  Liver, Lung, Prostate)
MultiPheno <- list(BreastPheno, ColonPheno, KidneyPheno,  LiverPheno, LungPheno, ProstatePheno)

MultiPheno <- lapply(MultiPheno, function(x){
  rownames(x) <- gsub('-','.',rownames(x))
  return(x)
})

for (i in 1:length(tissues)){
  
  temp <- intersect(colnames(MultiMat[[i]]), rownames(MultiPheno[[i]]))
  MultiMat[[i]] <- MultiMat[[i]][, temp]
  MultiPheno[[i]]   <- MultiPheno[[i]][temp,]
  
}

names(MultiMat) <- tissues
names(MultiPheno) <-  tissues

lapply(1:length(tissues), function(i){
  all(rownames(MultiPheno[[i]])== colnames(MultiMat[[i]]))
})


### check unique values
lapply(MultiMat, function(x){unique(unlist(x))})

### get common samples with IDList
MultiCnv <- lapply(1:length(MultiMat), function(i){
  Tempmat <- data.matrix(MultiMat[[i]])
  cnv <- Tempmat
  
  ### keep -2 and +2
  cnv[abs(cnv) < 2] <- 0
  
  return(cnv)
})

names(MultiCnv) <- tissues



### choose common genes and common sampels across different data type
tissues <- c('Breast','Colon', 'Kidney', 'Liver', 'Lung', 'Prostate')

CommonSample<- list()
CommonGenes <- list()
rna.combined<- list()
mut.combined<- list()
cnv.combined<- list()
pheno.combined <- list()

for(i in 1:length(tissues)){
  
  ### choose tumor samples
  cmsp <- Reduce(intersect, list(colnames(MultiDiv[[i]]), colnames(MultiMut[[i]]), colnames(MultiCnv[[i]])))
  cmgen <- Reduce(intersect, list(rownames(MultiDiv[[i]]), rownames(MultiMut[[i]]), rownames(MultiCnv[[i]])))
  
  CommonSample[[i]] <- cmsp
  CommonGenes[[i]] <- cmgen
  
  #### get phenotype
  
  temppheno <- MultiPheno[[i]][cmsp,]
  pheno.combined[[i]] <- temppheno
  
  #### build mat 
  
  rnamat <- MultiDiv[[i]][cmgen, cmsp]
  mutmat <-  MultiMut[[i]][cmgen, cmsp]
  cnvmat <- MultiCnv[[i]][cmgen, cmsp]

  rna.combined[[i]] <- rnamat
  mut.combined[[i]] <- mutmat
  cnv.combined[[i]] <- cnvmat
  
  
}

names(rna.combined) <- tissues
names(mut.combined) <- tissues
names(cnv.combined) <- tissues
names(pheno.combined) <- tissues


#####################################Construct Pairs from PathwayCommons#########

#####################################Reactome Pairs######################
#### read pairs 
Reactome <- read.csv(sprintf("%s/MECHANISM/PathwayCommons10.reactome.hgnc.txt.gz", DATA_DIR),
                     sep='\t', header = TRUE,stringsAsFactors = FALSE)
control_relation <-c("controls-expression-of", "controls-state-change-of")
PC_Reactome <- Reactome[Reactome$INTERACTION_TYPE %in% control_relation &Reactome$PATHWAY_NAMES !="" ,]


#### build table for weak connection and strong connection
StrongConnect <- PC_Reactome[PC_Reactome$INTERACTION_TYPE =='controls-expression-of', c(1:3)]
WeakConnect <- PC_Reactome[PC_Reactome$INTERACTION_TYPE=='controls-state-change-of', c(1:3)]

colnames(StrongConnect) <- c('source', 'interaction','target')
colnames(WeakConnect) <- c('source', 'interaction','target')

##### connect two genes with k=1
ConnectMatrix1 <- StrongConnect
ConnectMatrix1 <- cbind(ConnectMatrix1, 'S')[,c(1,3,4)]
colnames(ConnectMatrix1) <- c('source','target','interaction')
ConnectMatrixK1 <- ConnectMatrix1


#################### connect two genes with k<=2
#### weak-strong
ConnectMatrix21 <- inner_join(WeakConnect,StrongConnect, by=c('target' ='source'))
ConnectMatrix21 <- cbind(ConnectMatrix21, 'W-S')[,c(1,5,6)]
colnames(ConnectMatrix21) <- c('source','target','interaction')
#### strong-strong
ConnectMatrix22 <- inner_join(StrongConnect, StrongConnect, by=c('target'='source'))
ConnectMatrix22 <- cbind(ConnectMatrix22, 'S-S')[,c(1,5,6)]
colnames(ConnectMatrix22) <- c('source','target','interaction')
ConnectMatrix2 <- unique(rbind(ConnectMatrix21, ConnectMatrix22))


#### all connections with k<=2
ConnectMatrixK2<-unique(rbind(ConnectMatrix2, ConnectMatrixK1))

######################### connect two genes with k<=3
#### weak-weak-strong
Tempconnect31 <- inner_join(WeakConnect, WeakConnect, by=c('target' ='source'))
ConnectMatrix31 <- inner_join(Tempconnect31, StrongConnect, by=c('target.y' ='source'))
ConnectMatrix31 <- cbind(ConnectMatrix31,'W-W-S')[,c(1,7,8)]
colnames(ConnectMatrix31) <- c('source','target','interaction')
#### weak-strong-strong
Tempconnect32 <- inner_join(WeakConnect, StrongConnect, by=c('target' ='source'))
ConnectMatrix32 <- inner_join(Tempconnect32, StrongConnect, by=c('target.y' ='source'))
ConnectMatrix32 <- cbind(ConnectMatrix32,'W-S-S')[,c(1,7,8)]
colnames(ConnectMatrix32) <- c('source','target','interaction')
#### strong-strong-strong
Tempconnect33 <- inner_join(StrongConnect, StrongConnect, by=c('target' ='source'))
ConnectMatrix33 <- inner_join(Tempconnect32, StrongConnect, by=c('target.y' ='source'))
ConnectMatrix33 <- cbind(ConnectMatrix33,'S-S-S')[,c(1,7,8)]
colnames(ConnectMatrix33) <- c('source','target','interaction')
#### strong-weak-strong
Tempconnect34 <- inner_join(StrongConnect, WeakConnect, by=c('target' ='source'))
ConnectMatrix34 <- inner_join(Tempconnect34,StrongConnect, by=c('target.y' ='source'))
ConnectMatrix34 <- cbind(ConnectMatrix34,'S-W-S')[,c(1,7,8)]
colnames(ConnectMatrix34) <- c('source','target','interaction')

ConnectMatrix3 <- unique(rbind(ConnectMatrix31, ConnectMatrix32, ConnectMatrix33, ConnectMatrix34))

#### all connections with k<=3
ConnectMatrixK3<-unique(rbind(ConnectMatrix3, ConnectMatrixK2))
ConnectMatrix.list <- list(ConnectMatrixK1=ConnectMatrixK1, ConnectMatrixK2=ConnectMatrixK2, ConnectMatrixK3=ConnectMatrixK3)
names(ConnectMatrix.list) <- paste('k<=', c(1,2,3))

lapply(ConnectMatrix.list, dim)


#### remove pairs without available DNA/RNA informations
ConnectMatrixK3.combined <- lapply(1:length(mut.combined), function(i){
  mat_mut <- mut.combined[[i]]
  connectmatrix <- ConnectMatrix.list$`k<= 3`
  temp_id <- connectmatrix[,1] %in% rownames(mat_mut) & connectmatrix[,2] %in% rownames(mat_mut)
  temp <- connectmatrix[temp_id,]
  return(temp)
  
})

names(ConnectMatrixK3.combined) <- names(mut.combined)

Size_ConnectMatrixK3<-lapply(ConnectMatrixK3.combined, function(x){
  c(nrow(x), length(unique(x[,1])),length(unique(x[,2])) )
})

#### DNA aberration binary expresson
dna.combined <- lapply(1:length(mut.combined), function(i){
  mat_mut  <-abs(mut.combined[[i]])
  mat_cnv  <- abs(cnv.combined[[i]])
  temp_mat <- mat_mut + mat_cnv
  temp_mat[temp_mat >0] <- 1
  
  return(temp_mat)
  
})

names(dna.combined) <- names(rna.combined)



#### use chisq.test to compute independence test statistics 
#### independence test for all pairs within K<=3
Independence_K3.combined <- lapply(1:length(dna.combined), function(i){
  mat_dna <- dna.combined[[i]]
  mat_rna <- rna.combined[[i]]
  connectmatrix <- ConnectMatrixK3.combined[[i]]
  connectmatrix_sel <- connectmatrix[connectmatrix[,1] %in% rownames(mat_dna) & connectmatrix[,2] %in% rownames(mat_dna),]
  
  temp<-  apply(connectmatrix_sel, 1,function(x){
    
    x <- as.character(x)
    if (length(unique(mat_dna[x[1],])) >1 & length(unique(mat_rna[x[2],]))>1){
      t <- chisq.test(mat_dna[x[1],], mat_rna[x[2],])
      return(c(t$statistic, t$p.value, x[3]))
    }else{
      return(NULL)
    }
  })
  
  names(temp) <- paste(connectmatrix_sel[,1], connectmatrix_sel[,2], sep='_')
  temp <- do.call(rbind, temp)
  colnames(temp) <- c('chisq statistics','pvalue','interaction')
  return(temp)
  
})

names(Independence_K3.combined) <- names(dna.combined)


#### choose unique pairs with p value less than 0.05

Pairs_select_0.05 <- lapply(1:length(Independence_K3.combined), function(i){
  ind_sta<- Independence_K3.combined[[i]]
  pairs <- do.call(rbind, strsplit(rownames(ind_sta), split='_'))
  temp <- cbind(pairs, ind_sta)
  ### choose pvalue <=0.05
  temp_sel <- temp[as.numeric(temp[,4]) <= 0.05,]
  
  return(temp_sel)
})

lapply(Pairs_select_0.05, dim)
names(Pairs_select_0.05) <- names(Independence_K3.combined)

Pairs_select_0.05_unique <- lapply(Pairs_select_0.05, function(x){
  unique(x[,1:4])
})

lapply(Pairs_select_0.05_unique, dim)


### binary expression for pair level 
Motif_cmb.combined <- lapply(1:length(Pairs_select_0.05_unique), function(i){
  pairs_sel  <- Pairs_select_0.05_unique[[i]]
  mat_dna <- abs(dna.combined[[i]])
  mat_rna <- abs(rna.combined[[i]])
  Temp <- apply(pairs_sel, 1,function(x){
    x <- as.character(x)
    temp1 <- mat_dna[x[1],]
    temp2 <- mat_rna[x[2],]
    ### dna aberrant and rna divergence
    temp <- (temp1 ==1 & temp2==1)*1
    return(temp)
    
  })
  
  Temp <- t(Temp)
  return(Temp)
  
  
})
names(Motif_cmb.combined) <- tissues

### compute probability 
Divp_Motif_cmb.combined <- DivProb(Motif_cmb.combined)

### select pairs with div prob  >= 0.02
Pairs_select_0.05_unique_sel <- lapply(1:length(Divp_Motif_cmb.combined), function(i){
  
  sel_id <- names(Divp_Motif_cmb.combined[[i]])[Divp_Motif_cmb.combined[[i]] >=0.02]
  temp <- Pairs_select_0.05_unique[[i]][sel_id,]
  return(temp)
  
})
names(Pairs_select_0.05_unique_sel) <- tissues

Motif_cmb.combined_sel <- lapply(1:length(Motif_cmb.combined), function(i){
  sel_id <- which(Divp_Motif_cmb.combined[[i]] >=0.02)
  temp <- Motif_cmb.combined[[i]][sel_id,]
  return(temp)
})

names(Motif_cmb.combined_sel) <- tissues


#################################### binary expression for source#########################
### build source module on the pairs with pvalue <=0.05
source_module <- lapply(1:length(Pairs_select_0.05_unique), function(i){
  pairs <- Pairs_select_0.05_unique[[i]]
  source_gene <- unique(pairs[,1])
  
  Temp<-lapply(source_gene, function(x){
    temp <- rownames(pairs)[pairs[,1]==x]
    return(temp)
  })
  names(Temp) <- source_gene
  return(Temp)
  
})

names(source_module) <- names(Pairs_select_0.05_unique)

### binary expression for source motif
source_Motif_cmb.combined <- MotifExp(source_module, Motif_cmb.combined)
### divergent probability of source motif
divp_source_Motif_cmb.combined <- DivProb(source_Motif_cmb.combined)

#### select source genes with div prob >=0.02
source_module_sel <- lapply(1:length(source_module), function(i){
  sel_id <- names(source_module[[i]])[divp_source_Motif_cmb.combined[[i]] >=0.02]
  temp <- source_module[[i]][sel_id]
  return(temp)
})
names(source_module_sel) <- tissues
#### binary expression for source gene with div prob >= 0.02
source_Motif_cmb.combined_sel <- lapply(1:length(source_Motif_cmb.combined), function(i){
  sel_id <- rownames(source_Motif_cmb.combined[[i]])[divp_source_Motif_cmb.combined[[i]] >=0.02]
  temp <- source_Motif_cmb.combined[[i]][sel_id,]
  return(temp)
})
names(source_Motif_cmb.combined_sel) <- tissues




################################## binary expression for target######################


### build target module on the pairs with pvalue <=0.05
target_module <- lapply(1:length(Pairs_select_0.05_unique), function(i){
  pairs <- Pairs_select_0.05_unique[[i]]
  target_gene <- unique(pairs[,2])
  
  Temp<-lapply(target_gene, function(x){
    temp <- rownames(pairs)[pairs[,2]==x]
    return(temp)
  })
  names(Temp) <- target_gene
  return(Temp)
  
})

names(target_module) <- names(Pairs_select_0.05_unique)

### binary expression for target motif
target_Motif_cmb.combined <- MotifExp(target_module, Motif_cmb.combined)
### divergent probability of target motif
divp_target_Motif_cmb.combined <- DivProb(target_Motif_cmb.combined)
#### choose target module with div prob  >= 0.02
target_module_sel <- lapply(1:length(target_module), function(i){
  sel_id <- names(target_module[[i]])[divp_target_Motif_cmb.combined[[i]] >=0.02]
  temp <- target_module[[i]][sel_id]
  return(temp)
})
names(target_module_sel) <- tissues
#### binary expression for target gene with div prob >= 0.02
target_Motif_cmb.combined_sel <- lapply(1:length(target_Motif_cmb.combined), function(i){
  sel_id <- rownames(target_Motif_cmb.combined[[i]])[divp_target_Motif_cmb.combined[[i]] >=0.02]
  temp <- target_Motif_cmb.combined[[i]][sel_id,]
  return(temp)
})
names(target_Motif_cmb.combined_sel) <- tissues


###################################covering function ######################
pair_Result_Motif <- MinimalSigGene2(Motif_cmb.combined_sel, alpha = 0.00, J=1)
source_Result_Motif<- MinimalSigGene2(source_Motif_cmb.combined_sel, alpha = 0.00, J=1)
target_Result_Motif<- MinimalSigGene2(target_Motif_cmb.combined_sel, alpha = 0.00, J=3)



