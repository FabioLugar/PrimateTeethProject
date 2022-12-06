#HEADER-----
setwd("~/PrimateMolars/GlobalModels/")
library(mvtnorm)
library(PCMFit)
library(PCMBase)
library(PCMkappa)
library(data.table)
library(ape)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(plyr)
library(abind)
library(dplyr)
library(phytools)
library(psych)
library(evolqg)
library(reshape2)
library(treeplyr)
source('DefineParameterLimits_GLOBAL.R', local=FALSE)
doParallel::registerDoParallel(cores = 50)
options(PCMBase.Threshold.EV = 1e-8)
load("../sp_means_trees_updated.Rdata")

#Getting good Ps-----
tab4imcov<-dplyr::select(molarMeasures %>% subset(., M3!=0), Species, M2:M3) %>% na.omit
tab4imcov<-data.table(tab4imcov)[, if(.N > 30) .SD,Species]

Ps<-alply(unique(tab4imcov$Species), 1,
          function(i) {
            cov(tab4imcov[tab4imcov$Species==i,-1])
          })
names(Ps)<-unique(tab4imcov$Species)

Ps$pooled<-CalculateMatrix(lm(cbind(M2,M3)~Species,tab4imcov))

#Generating pseudoG-----
load("../Gbaboon2.Rdata")

Gs<-llply(1:10000,function(x){
  x<-rmvnorm(100, mean = Gmeans,sigma = Gcov)
  x<-x[,1:3]*x[,4:6]
  x<-x[,2:3]/x[,1]
  cov(x)
})
Ps$G<-apply(abind(Gs, along=3), 1:2, median)
dimnames(Ps$G)<-list(c("M2","M3"),c("M2","M3"))

#Matching tree to data-----
X<-means[,rownames(Ps$pooled)]
rownames(X)<-means$Species
#figuring out what to include
X<-X[means$Species!="Megaladapis_madagascariensis",]
X[X==0 &!is.na(X)]<-NA
X<-t(X)
X<-X[,!apply(X,2, function(x) all(is.na(x)))]

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
# now build models ----

modelNames<-PCMDefaultModelTypes()[c(2,4,5)]

modelNames<-setNames(modelNames,c("BM", "OUD", "OU"))
models<-llply(modelNames,PCM,k=2,modelTypes = modelNames,
              params = list(X0=c(1,1), 
                            Sigma_x=abind(D_x, along = 3)))
models$OU$Theta[] <- 1
models$OUD$Theta[]<- 1

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
                     "__Omitted_Sigmae_x"),
              paste0("OUkappa",
                     "__Global_X0",
                     "__NonNegative_kappa",
                     "__Diagonal_WithNonNegativeDiagonal_H",
                     "__Theta",
                     "__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Fixed_Sigma_x",
                     "__Omitted_Sigmae_x"))
modelNamesKappaP<-setNames(modelNamesKappa,c("BMkappaP", "OUkappaP", "OUkappaPdiag"))
modelNamesKappaG<-setNames(modelNamesKappa,c("BMkappaG", "OUkappaG", "OUkappaGdiag"))
modelNamesKappaD<-setNames(modelNamesKappa,c("BMkappaD", "OUkappaD", "OUkappaDdiag"))
models<-c(models,llply(modelNamesKappaP,PCM,k=2,modelTypes = modelNames,
              params = list(X0=c(1,1), 
                            kappa=1,
                            Sigma_x=abind(UpperTriFactor(Ps$pooled*0.5),
                                          along = 3))),
          llply(modelNamesKappaG,PCM,k=2,modelTypes = modelNames,
                params = list(X0=c(1,1), 
                              kappa=1,
                              Sigma_x=abind(UpperTriFactor(Ps$G),
                                            along = 3))),
          llply(modelNamesKappaD,PCM,k=2,modelTypes = modelNames,
                params = list(X0=c(1,1), 
                              kappa=1,
                              Sigma_x=abind(Ps$G%*%diag(2)%*%Ps$G*1000,
                                            along = 3))))
models$OUkappaP$Theta[]<-1
models$OUkappaG$Theta[]<-1
models$OUkappaD$Theta[]<-1
models$OUkappaPdiag$Theta[]<-1
models$OUkappaGdiag$Theta[]<-1
models$OUkappaDdiag$Theta[]<-1

SEsp<-matrix(0, nrow=2, ncol=PCMTreeNumTips(tree),dimnames = list(c("M2", "M3"),tree$tip.label))
M2errors<-dplyr::select(molarMeasures, Species, M2) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(M2), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(M2errors)<-M2errors$Species
M3errors<-dplyr::select(molarMeasures, Species, M3) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(M3), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(M3errors)<-M3errors$Species

for (i in tree$tip.label) SEsp[,i]<-unlist(c(M2errors[i,"se"],M3errors[i,"se"]))

SEsp[1,is.na(SEsp[1,])]<-weighted.mean(M2errors$sd,M2errors$n,na.rm = T)
SEsp[2,is.na(SEsp[2,])|SEsp[2,]==0]<-weighted.mean(M3errors$sd,M3errors$n,na.rm = T)

#Fitting full models --------

fitBM<-PCMFit(model = models$BM, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOU<-PCMFit(model = models$OU, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUdiag<-PCMFit(model = models$OUD, tree = tree, X = X,
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
fitBMkappaD<-PCMFit(model = models$BMkappaD, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaD<-PCMFit(model = models$OUkappaD, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaPdiag<-PCMFit(model = models$OUkappaPdiag, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaGdiag<-PCMFit(model = models$OUkappaGdiag, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)
fitOUkappaDdiag<-PCMFit(model = models$OUkappaDdiag, tree = tree, X = X,
              metaI = PCMBaseCpp::PCMInfoCpp,
              SE = SEsp, doParallel =T)

save.image("fitGlobal.Rdata")
