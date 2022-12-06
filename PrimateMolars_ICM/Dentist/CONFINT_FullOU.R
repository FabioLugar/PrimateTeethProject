#HEADER------
setwd("~/PrimateMolars/Dentist/")
library(PCMBase)
library(PCMkappa)
library(PCMFit)
library(dentist)
library(ggplot2)
library(ggbeeswarm)
library(reshape)
library(cowplot)
library(foreach)
library(plyr)
library(abind)
library(dplyr)
library(expm)
library(psych)
options(PCMBase.Threshold.EV = 1e-8)

load("../GlobalModels/fitGlobal.Rdata")

bestModel<-RetrieveBestModel(fitOU)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0", "H1", "H2", "H3","theta1","theta2","sigma1","sigma2", "sigma3")
best_neglnL<- c(-fitOU$logLikOptim)
# PCMLik(X,tree, bestModel, SEsp, metaI=PCMBaseCpp::PCMInfoCpp)

wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.01,-0.5,0.01,0.8,0.8,0.01,0.01,0.01)
maxPars[]<-c(1.2,1.6,0.5 , 0.5,0.5 ,1.2,1.6,0.5 ,0.5 ,0.5 )

dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

dented_results$all_ranges[,c("sigma1","sigma2", "sigma3")]

H<-fitOU$modelOptim$H

sigmaG<-PCMApplyTransformation(fitOUkappaG$modelOptim)$Sigma_x[,,1] %>% tcrossprod

hxs<-dented_results$results %>% subset(., neglnL<(best_neglnL+8)) %>% unique %>%
  ddply(., 1, function(x) {
    s<-h<-matrix(0,2,2)
    s[upper.tri(s,diag = T)]<-unlist(x[,c("sigma1","sigma2","sigma3")])
    s<-tcrossprod(s)
    H[,,1][upper.tri(s,diag = T)]<-unlist(x[,c("H1","H2","H3")])
    H<-PCMApplyTransformation(H)
    data.frame(sigma=tr(s),H=tr(H[,,1]), deltaL=abs(best_neglnL-x[,"neglnL"]))
  })

ggplot(hxs, aes(sigma, H))+
  geom_point()+
  geom_smooth(method="lm")+
  xlab(expression(tr(Sigma)))+
  ylab(expression(tr(H)))