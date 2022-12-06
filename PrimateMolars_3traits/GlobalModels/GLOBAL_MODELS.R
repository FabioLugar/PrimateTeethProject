#HEADER-----
setwd("~/PrimateMolars_3traits/GlobalModels/")
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
Ps$G<-apply(abind(Gs, along=3), 1:2, mean)
dimnames(Ps$G)<-list(c("m1a","m2a","m3a"),c("m1a","m2a","m3a"))

#Matching tree to data-----
X<-log(means[,rownames(Ps$pooled)])
rownames(X)<-means$Species
#figuring out what to include

X<-X[!apply(select(means,M2,M3),1, function(x) all(is.na(x))) & means$Species!="Megaladapis_madagascariensis",]
X[X==-Inf &!is.na(X)]<-NA
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

modelNames<-PCMDefaultModelTypes()[c(2,5)]
modelNames<-setNames(modelNames,c("BM", "OU"))
models<-llply(modelNames,PCM,k=3,modelTypes = modelNames,
              params = list(X0=c(1,1,1), 
                            Sigma_x=abind(D_x, along = 3)))
models$OU$Theta[] <- 1

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
              params = list(X0=c(1,1,1), 
                            kappa=1,
                            Sigma_x=abind(UpperTriFactor(Ps$pooled*0.5),
                                          along = 3))),
          llply(modelNamesKappaG,PCM,k=3,modelTypes = modelNames,
                params = list(X0=c(1,1,1), 
                              kappa=1,
                              Sigma_x=abind(UpperTriFactor(Ps$G),
                                            along = 3))))
models$OUkappaP$Theta[]<-1
models$OUkappaG$Theta[]<-1

SEsp<-matrix(0, nrow=3, ncol=PCMTreeNumTips(tree),dimnames = list(c("m1a","m2a", "m3a"),tree$tip.label))
m1aerrors<-dplyr::select(molarMeasures, Species, m1a) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m1a)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m1aerrors)<-m1aerrors$Species
m2aerrors<-dplyr::select(molarMeasures, Species, m2a) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m2a)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m2aerrors)<-m2aerrors$Species
m3aerrors<-dplyr::select(molarMeasures, Species, m3a) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m3a)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m3aerrors)<-m3aerrors$Species

for (i in tree$tip.label) SEsp[,i]<-unlist(c(m1aerrors[i,"se"],m2aerrors[i,"se"],m3aerrors[i,"se"]))

SEsp[1,is.na(SEsp[1,])]<-weighted.mean(m1aerrors$sd,m1aerrors$n,na.rm = T)
SEsp[2,is.na(SEsp[2,])]<-weighted.mean(m2aerrors$sd,m2aerrors$n,na.rm = T)
SEsp[3,is.na(SEsp[3,])|SEsp[3,]==0]<-weighted.mean(m3aerrors$sd,m3aerrors$n,na.rm = T)

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
