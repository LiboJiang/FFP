


pheno1 <- read.csv("Data/Daily Transpiration2.csv")

pheno_n <- pheno1[,1]
pheno2 <- matrix(as.numeric(as.matrix(pheno1[,-1])),nrow=245)


marker <- read.csv("Data/Table_S2-SNP-To Libo.csv")

marker_info <- marker[,1:3]
marker_m <- marker[,-c(1:3)]


ncom <- colnames(marker_m)

coid <- names(which(table(c(ncom,as.character(pheno_n)))>1))

pii <- c()
mii <- c()
for(i in 1:length(coid)){
  tmp1 <- which(pheno_n==coid[i])
  tmp2 <- which(ncom==coid[i])
  if(length(tmp2)>0){
    tmp3 <- rep(tmp2,length(tmp1))
  }else{
    next;
  }
  pii <- c(pii,tmp1);mii <- c(mii,tmp3)
}

npheno <- pheno2[pii,]
mmarker <- marker_m[,mii]

deli <- c()
for(i in 1:dim(mmarker)[1]){
  dl <- table(as.numeric(mmarker[i,]))
  if(length(dl)<=1){
    deli <- c(deli,i)
  }
}

mmarker1 <- mmarker[-deli,]
mmarker_info <- marker_info[-deli,]

dat <- list(phenotype=npheno,genotype=mmarker1,info=mmarker_info)
save(dat,file="dat.RData")






