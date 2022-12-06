setwd("~/PrimateMolars/")
library(ape)
library(phytools)
library(PCMFit)
library(PCMBase)
library(data.table)
library(foreach)
library(doParallel)
library(plyr)
library(dplyr)
library(geomorph)
library(ggrepel)
library(ggplot2)
library(cowplot)
library(ggtree)
library(RColorBrewer)
library(psych)
library(geiger)
library(mvtnorm)
library(ggstance)
library(evolqg)
source('DefineParameterLimits_shifts.R', local=FALSE)
doParallel::registerDoParallel(cores = 50)

load("fitMIXED.RData")
load("fitMixed5spBEST.RData")
treeFixedShifts <- PCMTree(tree) 
PCMTreeSetPartRegimes(
  treeFixedShifts, 
  part.regime = c(`481` = 1, `483` = 2,`651`= 3, `505`=1, `540`=3), 
  setPartition = TRUE)

PCMTreePlot(treeFixedShifts)

sp_regime<-data.frame(regime=treeFixedShifts$edge.part,treeFixedShifts$edge) %>%
  subset(., X2<as.numeric(treeFixedShifts$node.label[1]))
sp_regime$X2<-tree$tip.label[sp_regime$X2]
sp_regime<-sp_regime[,-2]
levels(sp_regime$regime)<-treeFixedShifts$part.regime[levels(sp_regime$regime)]
sp_regime$regime<-factor(sp_regime$regime,levels = 1:length(unique(sp_regime$regime)))

df2plot<-data.frame(t(X)[sp_regime$X2,],sp_regime)
df2plot[df2plot$regime==3 & is.na(df2plot$M3),"M3"]<-0

ggplot(df2plot, aes(M2, M3))+
  geom_point(aes(x=M2,y=M3, color=regime))+
  geom_smooth(aes(color=regime), subset(df2plot,M3>0),method="lm",se=F)+
  theme_bw()+
  xlab("M2/M1")+
  ylab("M3/M1")

mixedMod<-MixedGaussian(k=2,modelTypes=MGPMDefaultModelTypes()[5], 
                        mapping=c(`1`=1,`2`=1,`3`=1), 
                        X0=RetrieveBestModel(modelFits$OU)$X0,
                        Sigmae_x = structure(0, 
                                             class = c("MatrixParameter", "_Omitted"), 
                                             description = "upper triangular factor of the non-phylogenetic variance-covariance"))

bestModel<-RetrieveBestFitScore(fitMIXED[[9]])
mixedMod$`1`$H<-mixedMod$`2`$H<-mixedMod$`3`$H<-bestModel$inferredModel$`2`$H
mixedMod$`1`$Theta<-mixedMod$`2`$Theta<-mixedMod$`2`$Theta<-bestModel$inferredModel$`2`$Theta
mixedMod$`1`$Sigma_x<-mixedMod$`2`$Sigma_x<-mixedMod$`2`$Sigma_x<-bestModel$inferredModel$`2`$Sigma_x

testModel<-mixedMod
partrace<-NULL
for (i in 1:iter){
  modelFit<-PCMFit(model = testModel, tree = tree, X = X,
                   metaI = PCMBaseCpp::PCMInfoCpp,
                   numRunifInitVecParams = 1000,
                   numGuessInitVecParams = 100,
                   SE = SEsp, doParallel =T, numCallsOptim=50)
  
  testModel<-RetrieveBestModel(modelFit)
  logLiktrace[i]<-modelFit$logLikOptim[1]
  partrace<-rbind(partrace,modelFit$Optim$par)
  
  if(i %in% seq(100,iter,by=20)) {
    
    p1<-
      data.frame(i=1:i, logLik=logLiktrace[1:i]) %>%
      ggplot(., aes(i, logLik))+
      geom_path()+
      geom_point()+
      xlab("iter")
    
    p2<-
      data.frame(i=1:(i-1), dif=0.0001+logLiktrace[2:i]-logLiktrace[1:(i-1)]) %>%
      ggplot(., aes(i, dif))+
      geom_path()+
      geom_point()+
      scale_y_log10()+
      xlab("iter")
    
    p3<-partrace %>% melt(.) %>%
      ggplot(., aes(Var1, value))+
      facet_wrap(Var2~.,scales = "free_y")+
      geom_path()+
      geom_point()+
      xlab("iter")
    
    
    ggsave(paste0("logLikprofile_modelfitMixedFix_i",i,".pdf"),plot_grid(p1,p2,ncol = 1,align = "hv"),width = 7,height = 7)
    ggsave(paste0("paramsprofile_modelfitMixedFix_i",i,".pdf"),p3,width = 10,height = 7)
    save.image(paste0("modelfitMixedFix_i",i,".Rdata"))
    
    if(i!=100) {
      do.call(file.remove, list(list.files(pattern = paste0("modelfitMixedFix_i",i-20))))
    }
  }
}

tracelist$logLik$MixedFix<-logLiktrace
tracelist$param$MixedFix<-partrace
modelFits$MixedFix<-modelFit
do.call(file.remove, list(list.files(pattern = "modelfitMixedFix_i")))

save.image("fitGlobalBM_OU_MixedFixed.Rdata")
# 
# 
# names(fitMIXED)<-paste("Mixed",1:length(fitMIXED))
# load("fitMixed5spBEST.RData")
# fitMIXED$Best<-fitMONKEYS5sp
# 
# rbind(ldply(modelFits,function(mod) data.frame(loglik=logLik(mod), AIC=AIC(mod),BIC=BIC(mod))),
#       data.frame(.id="Fixed",loglik=logLik(testModel), AIC=AIC(testModel),BIC=BIC(testModel)))
# 
# ldply(fitMIXED,function(x){
#   mod<-RetrieveBestFitScore(x)$inferredModel
#   data.frame(loglik=logLik(mod), AIC=AIC(mod),BIC=BIC(mod))}),
# 
