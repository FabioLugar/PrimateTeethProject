#HEADER-----
setwd("~/PrimateMolars_3traits/GlobalModels_PCA/")
library(mvtnorm)
library(PCMFit)
library(PCMBase)
library(PCMkappa)
library(data.table)
library(ape)
library(ggtree)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(plyr)
library(abind)
library(dplyr)
library(evolqg)
library(phytools)
library(psych)
library(reshape2)
library(treeplyr)
source('DefineParameterLimits_GLOBAL.R', local=FALSE)
doParallel::registerDoParallel(cores = 50)
options(PCMBase.Threshold.EV = 1e-8)
load("../sp_means_trees_updated.Rdata")
means<-read.csv("imputed_for_pca.csv")
rownames(means)<-means$X
means<-means[,-1]
#Getting good Ps-----
tab4imcov<-dplyr::select(molarMeasures %>% subset(., M3!=0), Species, m1a:m3a) %>% na.omit
tab4imcov[,-1]<-log(tab4imcov[,-1])
tab4imcov<-data.table(tab4imcov)[, if(.N > 30) .SD,Species]

Ps<-alply(unique(tab4imcov$Species), 1,
          function(i) {
            cov(tab4imcov[tab4imcov$Species==i,-1])
          })
names(Ps)<-unique(tab4imcov$Species)

Ps$pooled<-CalculateMatrix(lm(cbind(m1a,m2a,m3a)~Species,tab4imcov))

#Generating pseudoG-----
load("../Gbaboon2.Rdata")

Gs<-llply(1:10000,function(x){
  x<-rmvnorm(100, mean = Gmeans,sigma = Gcov)
  x<-x[,1:3]*x[,4:6]
  cov(log(x))
})
Ps$G<-apply(abind(Gs, along=3), 1:2, median)
dimnames(Ps$G)<-list(c("m1a","m2a","m3a"),c("m1a","m2a","m3a"))

#Matching tree to data-----
X<-means[,rownames(Ps$pooled)]
#figuring out what to include
X<-X[rownames(means)!="Megaladapis_madagascariensis",]
# X[X==-Inf &!is.na(X)]<-NA
pca<-prcomp(X)
X<-pca$x

X<-t(X)
# X<-X[,!apply(X,2, function(x) all(is.na(x)))]

tree <- drop.tip(tree,tree$tip.label[!tree$tip.label %in% colnames(X)])
tree <- PCMTree(tree)
X<-X[,tree$tip.label]


X_na<-na.omit(t(X))
tree_na<-drop.tip(tree,
                  tree$tip.label[!tree$tip.label %in% rownames(X_na)]) %>%
  multi2di
pics<-apply(X_na,2,pic,phy=tree_na)
D<-crossprod(pics)/Nnode(tree_na)
D_x<-UpperTriFactor(D)


# error ----

SEsp<-array(NA,dim = c(3,3,PCMTreeNumTips(tree)), dimnames = list(rownames(X),rownames(X),colnames(X)))

for(i in names(Ps)){
  if(i %in% dimnames(SEsp)[[3]]){
    SEsp[,,i]<- UpperTriFactor((t(pca$rotation) %*% Ps[[i]] %*% pca$rotation)/table(molarMeasures$Species)[i]) 
  }
}

for(i in 1:dim(SEsp)[3]){
  if(is.na(sum(SEsp[,,i]))){
    SEsp[,,i] <- UpperTriFactor(t(pca$rotation) %*% Ps$pooled %*% pca$rotation) 
  }
}


# now build models ----

modelNames<-PCMDefaultModelTypes()[c(2,5)]
modelNames<-setNames(modelNames,c("BM", "OU"))
models<-llply(modelNames,PCM,k=3,modelTypes = modelNames,
              params = list(X0=c(0,0,0), 
                            Sigma_x=abind(D_x, along = 3)))
models$OU$Theta[] <- 0

modelNamesKappa<-c(paste0("BMkappa",
                     "__Global_X0",
                     "__NonNegative_kappa",
                     "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Fixed_Sigma_x",
                     "__Omitted_Sigmae_x"),
              paste0("OUkappa",
                     "__Global_X0",
                     "__NonNegative_kappa",
                     "__Schur_UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Transformable_H",
                     "__Theta",
                     "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Fixed_Sigma_x",
                     "__Omitted_Sigmae_x"))
modelNamesKappaP<-setNames(modelNamesKappa,c("BMkappaP", "OUkappaP"))
modelNamesKappaG<-setNames(modelNamesKappa,c("BMkappaG", "OUkappaG"))
models<-c(models,llply(modelNamesKappaP,PCM,k=3,modelTypes = modelNames,
              params = list(X0=c(0,0,0), 
                            kappa=1,
                            Sigma_x=abind(UpperTriFactor(t(pca$rotation) %*% (Ps$pooled*0.5) %*% pca$rotation),
                                          along = 3))),
          llply(modelNamesKappaG,PCM,k=3,modelTypes = modelNames,
                params = list(X0=c(0,0,0), 
                              kappa=1,
                              Sigma_x=abind(UpperTriFactor(t(pca$rotation) %*% Ps$G %*% pca$rotation),
                                            along = 3))))
models$OUkappaP$Theta[]<-0
models$OUkappaG$Theta[]<-0



#Fitting full models --------

fitBM<-PCMFit(model = models$BM, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)

fitOU<-PCMFit(model = models$OU, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitBMkappaP<-PCMFit(model = models$BMkappaP, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaP<-PCMFit(model = models$OUkappaP, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitBMkappaG<-PCMFit(model = models$BMkappaG, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaG<-PCMFit(model = models$OUkappaG, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)

save.image("fitGlobal.Rdata")
