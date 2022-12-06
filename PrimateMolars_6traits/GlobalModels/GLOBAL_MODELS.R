#HEADER-----
setwd("~/PrimateMolars_6traits/GlobalModels/")
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
tab4imcov<-dplyr::select(molarMeasures %>% subset(., M3!=0), Species, m1bl:m3md) %>% na.omit
tab4imcov[,-1]<-log(tab4imcov[,-1])
tab4imcov<-data.table(tab4imcov)[, if(.N > 30) .SD,Species]

Ps<-alply(unique(tab4imcov$Species), 1,
          function(i) {
            cov(tab4imcov[tab4imcov$Species==i,-1])
          })
names(Ps)<-unique(tab4imcov$Species)

Ps$pooled<-CalculateMatrix(lm(cbind(m1bl,m1md,m2bl,m2md,m3bl,m3md)~Species,tab4imcov))

#Generating pseudoG-----
load("../Gbaboon2.Rdata")

Gs<-llply(1:10000,function(x){
  x<-rmvnorm(100, mean = Gmeans,sigma = Gcov)
  cov(log(x))
})
Ps$G<-apply(abind(Gs, along=3), 1:2, mean)
dimnames(Ps$G)<-list(c("m1bl","m1md","m2bl","m2md","m3bl","m3md"),
                     c("m1bl","m1md","m2bl","m2md","m3bl","m3md"))

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
models<-llply(modelNames,PCM,k=6,modelTypes = modelNames,
              params = list(X0=c(1,1,1,1,1,1), 
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
models<-c(models,llply(modelNamesKappaP,PCM,k=6,modelTypes = modelNames,
              params = list(X0=c(1,1,1,1,1,1), 
                            kappa=1,
                            Sigma_x=abind(UpperTriFactor(Ps$pooled*0.5),
                                          along = 3))),
          llply(modelNamesKappaG,PCM,k=6,modelTypes = modelNames,
                params = list(X0=c(1,1,1,1,1,1), 
                              kappa=1,
                              Sigma_x=abind(UpperTriFactor(Ps$G),
                                            along = 3))))
models$OUkappaP$Theta[]<-1
models$OUkappaG$Theta[]<-1

SEsp<-matrix(0, nrow=6, ncol=PCMTreeNumTips(tree),dimnames = list(c("m1bl","m1md","m2bl","m2md","m3bl","m3md"),tree$tip.label))
m1blerrors<-dplyr::select(molarMeasures, Species, m1bl) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m1bl)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m1blerrors)<-m1blerrors$Species
m1mderrors<-dplyr::select(molarMeasures, Species, m1md) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m1md)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m1mderrors)<-m1mderrors$Species

m2blerrors<-dplyr::select(molarMeasures, Species, m2bl) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m2bl)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m2blerrors)<-m2blerrors$Species
m2mderrors<-dplyr::select(molarMeasures, Species, m2md) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m2md)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m2mderrors)<-m2mderrors$Species

m3blerrors<-dplyr::select(molarMeasures, Species, m3bl) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m3bl)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m3blerrors)<-m3blerrors$Species
m3mderrors<-dplyr::select(molarMeasures, Species, m3md) %>% 
  group_by(., Species) %>% na.omit %>%
  summarise(., sd=sd(log(m3md)), n=n()) %>% 
  mutate(., se=sd/sqrt(n)) %>% as.data.frame
rownames(m3mderrors)<-m3mderrors$Species

for (i in tree$tip.label) SEsp[,i]<-unlist(c(m1blerrors[i,"se"],m1mderrors[i,"se"],
                                             m2blerrors[i,"se"],m2mderrors[i,"se"],
                                             m3blerrors[i,"se"],m3mderrors[i,"se"]))

SEsp[1,is.na(SEsp[1,])]<-weighted.mean(m1blerrors$sd,m1blerrors$n,na.rm = T)
SEsp[2,is.na(SEsp[2,])]<-weighted.mean(m1mderrors$sd,m1mderrors$n,na.rm = T)
SEsp[3,is.na(SEsp[3,])]<-weighted.mean(m2blerrors$sd,m2blerrors$n,na.rm = T)
SEsp[4,is.na(SEsp[4,])]<-weighted.mean(m2mderrors$sd,m2mderrors$n,na.rm = T)
SEsp[5,is.na(SEsp[5,])|SEsp[5,]==0]<-weighted.mean(m3blerrors$sd,m3blerrors$n,na.rm = T)
SEsp[6,is.na(SEsp[6,])|SEsp[6,]==0]<-weighted.mean(m3mderrors$sd,m3mderrors$n,na.rm = T)

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
