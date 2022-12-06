setwd("~/PrimateMolars/EXPLORING_PCM/")
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
library(reshape)
library(expm)
options(PCMBase.Threshold.EV = 1e-8)

doParallel::registerDoParallel(cores = 50)
TipLabelTable <- function(model, ...) {
  tree <- PCMTree(attr(model, "tree"))
  data.table(
    node = sapply(sapply(PCMRegimes(tree), 
                         PCMTreeGetTipsInRegime, 
                         tree = tree), sample, size = 1),
    part.model = 
      paste0(" ", PCMRegimes(tree), ".", 
             LETTERS[PCMMapModelTypesToRegimes(model, tree)], " "), 
    ...)
}


#Global model----
load("~/PrimateMolars/GlobalModels/fitGlobal.Rdata")

ggplot(subset(molarMeasures, Species!="Megaladapis_madagascariensis"), aes(M2,M3))+
  scale_fill_brewer(name="",palette = "Set3")+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  geom_point(alpha=1)+
  theme_minimal()+
  coord_cartesian()+
  xlab("M2/M1")+
  ylab("M3/M1")
  # coord_fixed()


bestModels<-list(BM=RetrieveBestModel(fitBM),
                 OU=RetrieveBestModel(fitOU),
                 OUD=RetrieveBestModel(fitOUdiag),
                 BMkappaP=RetrieveBestModel(fitBMkappaP),
                 OUkappaP=RetrieveBestModel(fitOUkappaP),
                 BMkappaG=RetrieveBestModel(fitBMkappaG),
                 OUkappaG=RetrieveBestModel(fitOUkappaG),
                 BMkappaD=RetrieveBestModel(fitBMkappaD),
                 OUkappaD=RetrieveBestModel(fitOUkappaD),
                 OUkappaPdiag=RetrieveBestModel(fitOUkappaPdiag),
                 OUkappaGdiag=RetrieveBestModel(fitOUkappaGdiag),
                 OUkappaDdiag=RetrieveBestModel(fitOUkappaDdiag))

modelCompare<-
  list(BM=fitBM,
       OU=fitOU,
       OUdiag=fitOUdiag,
       BMkappaP=fitBMkappaP,
       OUkappaP=fitOUkappaP,
       BMkappaG=fitBMkappaG,
       OUkappaG=fitOUkappaG,
       OUkappaGdiag=fitOUkappaGdiag) %>%
  ldply(.,function(x) data.frame(logLik=x$logLikOptim,
                                 n=dim(x$X)[2],
                                 k=PCMParamCount(x$modelOptim)))  %>%
  mutate(
    .,
    BIC = (log(n) * k - 2 * logLik)
  )

arrange(modelCompare, BIC)

bestModel<-bestModels$OUkappaG

stationary <-PCMkappa:::StationaryVariance(bestModel$H[,,1],bestModel$Sigma_x[,,])
nonPhylovar<-var(na.omit(t(X)))

rotH<-t(eigen(stationary)$vectors) %*% PCMApplyTransformation(bestModel)$H[,,1] %*% eigen(stationary)$vectors
log(2)/diag(rotH)
max(nodeHeights(tree))

log(2)/diag(PCMApplyTransformation(bestModel)$H[,,1])

statElli<-ellipse::ellipse(stationary, centre=bestModel$Theta[])
nphyElli<-ellipse::ellipse(nonPhylovar,centre=colMeans(t(X),na.rm = T))
ancGElli<-ellipse::ellipse(Ps$pooled*0.5,centre=bestModel$X0[])
baboonGElli<-ellipse::ellipse(Ps$G,centre=bestModel$X0[])
rateElli<-ellipse::ellipse(tcrossprod(PCMApplyTransformation(bestModel)$Sigma_x[,,]),centre=bestModel$X0[])


H<-PCMApplyTransformation(bestModel)$H[,,1]
Sigma<-tcrossprod(PCMApplyTransformation(bestModel)$Sigma_x[,,])
selectElli<-ellipse::ellipse(
  # solve(solve(sqrtm(Sigma))%*%H%*%solve(sqrtm(Sigma))),centre=bestModel$Theta[])
  solve(sqrtm(H))%*%Sigma%*%solve(sqrtm(H)) - Ps$pooled,centre=bestModel$Theta[])

ellis<-data.frame(rbind(nphyElli,statElli,rateElli,ancGElli,baboonGElli,selectElli), type=rep(c("Disparity","Evolutionary (stationary)","Sigma (stochastic)","G (P*0.5)", "MC-G","H (deterministic)"), each=100)) %>%
  mutate(., type=factor(type,unique(type)))

allEllis.plot<-
  ggplot(subset(means, Species!="Megaladapis_madagascariensis"), aes(M2,M3))+
  geom_point(color="gray",data=molarMeasures, size=3)+
  geom_point(size=1)+
  # coord_fixed()+
  geom_polygon(aes(M2,M3, fill=type),color="black",ellis, alpha=0.5)+
  scale_fill_brewer(palette = "Set2")
allEllis.plot

someEllis.plot<-
  ggplot(subset(means, Species!="Megaladapis_madagascariensis"), aes(M2,M3))+
  # coord_fixed()+
  geom_polygon(aes(M2,M3, fill=type),color="black",
               subset(ellis, type %in% c("Evolutionary (stationary)","Sigma (stochastic)","H (deterministic)")) %>% 
                 mutate(.,type=factor(type, levels=levels(type)[c(2,6,3)])), alpha=0.7)+
  scale_fill_brewer(name="",palette = "Set3")+
  geom_point(shape=21,alpha=0.2, fill="lightgray")+
  theme_minimal()+
  coord_cartesian()+
  xlab("M2/M1")+
  ylab("M3/M1")+
  coord_fixed()
someEllis.plot
ggsave("allEllis.pdf", allEllis.plot,width = 8,height = 7)
ggsave("someEllis.pdf",someEllis.plot,width = 6,height = 6)

#Mixed model----
load("~/PrimateMolars/MixedModels/BIC_results/fitMIXED.RData")
load("~/PrimateMolars/MixedModels_ShiftFixed/fitMIXED_fixed.RData")
names(fitMIXED)<-paste("Mixed",1:length(fitMIXED))

modelCompare<-
  list(BM=fitBM,
       OU=fitOU,
       BMkappaP=fitBMkappaP,
       OUkappaP=fitOUkappaP,
       BMkappaG=fitBMkappaG,
       OUkappaG=fitOUkappaG,
       MixedBMBMBM=fitBMBMBM,
       MixedOUBMOU=fitOUBMOU,
       MixedOUOUOU=fitOUOUOU) %>%
  ldply(.,function(x) data.frame(logLik=x$logLikOptim,
                                 n=dim(x$X)[2],
                                 k=PCMParamCount(x$modelOptim)))  %>%
  mutate(
    .,
    BIC = (log(n) * k - 2 * logLik)
  )

names(fitMIXED)<-paste("fitMixed_BIC",1:10)

modelCompare<-
  rbind(modelCompare,
        ldply(fitMIXED,function(x){
          mod<-RetrieveBestFitScore(x)$inferredModel
          data.frame(logLik=x$tableFitsRR$logLik[2],
                     n=x$tableFitsRR$nobs[2],
                     k=PCMParamCount(mod))
        }) %>% mutate(., BIC = (log(n) * k - 2 * logLik)))
arrange(modelCompare[,c(1,4,2,5)], BIC)
library(xtable)
xtable(modelCompare[,c(1,4,2,5)])

write.csv(modelCompare,"modelCompare.csv")

print(PCMTable(fitOUOUOU$modelOptim), xtable = TRUE)

# bestFitMixed<-RetrieveBestFitScore(fitMIXED$`fitMixed 8`)
# 
# startree<-stree(length(tree$tip.label), "star")
# startree$edge.length<-rep(max(nodeHeights(tree)), times=length(tree$tip.label))
# X0<-bestFitMixed$inferredModel$X0
# 
# SimTr<-ldply(bestFitMixed$inferredModel[-1], function(x){
#   model<-PCM(sub("Omitted_X0","Global_X0",class(x)[1]), k = 2)
#   model$X0<-X0
#   for(i in names(x)) model[i]<-x[i]
#   ldply(1:100, function(i) {
#     XX<-PCMSim(startree, model, X0=X0)
#     tr(var(t(XX)))
#   },.parallel = T)
# })
# 
# TipLabelTable <- function(model, ...) {
#   tree <- PCMTree(attr(model, "tree"))
#   data.table(
#     node = sapply(sapply(PCMRegimes(tree), 
#                          PCMTreeGetTipsInRegime, 
#                          tree = tree), sample, size = 1),
#     part.model = 
#       paste0(" ", PCMRegimes(tree), ".", 
#              LETTERS[PCMMapModelTypesToRegimes(model, tree)], " "), 
#     ...)
# }
# 
# attr(bestFitMixed$inferredModel, "tree")$tip.label<-tree$tip.label
# 
# tree.plot <- 
#   PCMTreePlot(attr(bestFitMixed$inferredModel, "tree")) %<+% 
#   TipLabelTable(bestFitMixed$inferredModel, offset = 5) + 
#   # geom_tiplab2(size = 2) + 
#   geom_tiplab(aes(label = part.model), offset = 25) + 
#   geom_tiplab(size=2)+
#   # geom_nodelab(size = 2, color = "black") + 
#   theme(legend.position = "none")+
#   xlim_tree(100)
# 
# tr.plot<-
#   ggplot(SimTr %>% mutate(.,.id=factor(.id, 1:10)), aes(.id, V1))+
#   geom_boxplot(aes(fill=.id), show.legend = F)+
#   # scale_y_log10()+
#   xlab("model")+
#   ylab("Trace")+
#   theme_minimal()
# 
# trait.plot<-
#   PCMPlotTraitData2D(attr(bestFitMixed$inferredModel, "X"),
#                    attr(bestFitMixed$inferredModel, "tree"),
#                    alphaPoints = 1,scaleSizeWithTime = F,
#                    numTimeFacets = 1)+
#   xlab("M2/M1")+
#   ylab("M3/M1")
# 
# slope_ints<-llply(bestFitMixed$inferredModel[-1], function(x){
#   if(class(x)[2]=="OU") v<-PCMkappa:::StationaryVariance(x$H[,,], x$Sigma_x[,,])
#   if(class(x)[2]=="BM") v<-tcrossprod(x$Sigma_x[,,])
#   eigen(v)$vector[,1]
# })
# 
# modelTree<-bestFitMixed$tree
# sp_regime<-data.frame(regime=modelTree$edge.part,modelTree$edge) %>%
#   subset(., X2<as.numeric(modelTree$node.label[1]))
# sp_regime$X2<-tree$tip.label[sp_regime$X2]
# sp_regime<-sp_regime[,-2]
# levels(sp_regime$regime)<-modelTree$part.regime[levels(sp_regime$regime)]
# sp_regime$regime<-factor(sp_regime$regime,levels = 1:length(unique(sp_regime$regime)))
# 
# df2plot<-data.frame(t(X)[sp_regime$X2,],sp_regime)
# df2plot[df2plot$regime==5 & is.na(df2plot$M3),"M3"]<-0
# 
# means<-group_by(df2plot, regime) %>% summarise_each(., function(x) mean(x,na.rm=T))
# 
# ggplot(df2plot, aes(M2, M3))+
#   geom_point(aes(x=M2,y=M3, color=regime))+
#   # geom_smooth(aes(color=regime),method="lm",se=F)+
#   theme_bw()+
#   xlab("M2/M1")+
#   ylab("M3/M1")
# 
# ggsave("Mixed_model_vars.pdf",plot_grid(tree.plot,trait.plot,tr.plot,ncol = 2),width = 10,height = 10)
# 
# 
# #########
# 
# names(fitMIXED)<-paste("Mixed",1:length(fitMIXED))
# rbind(ldply(fitMIXED,function(x){
#   mod<-RetrieveBestFitScore(x)$inferredModel
#   data.frame(loglik=logLik(mod), AIC=AIC(mod),BIC=BIC(mod))
# }))
# # modelFits
# # PCMTable()
# 
# tree.plots<-
#   llply(fitMIXED[1:sum(!laply(fitMIXED,is.null))], function(x){
#   bestFit <- RetrieveBestFitScore(x)
#   modelTree<-bestFit$tree
#   
#   modelTree$tip.label<-tree$tip.label
#   PCMTreePlot(modelTree)%<+% 
#     TipLabelTable(bestFit$inferredModel, offset = 5) + 
#     # geom_tiplab2(size = 2) + 
#     geom_tiplab(aes(label = part.model), offset = 25) + 
#     geom_tiplab(aes(label=sub("_"," ",label)), size=2, hjust=-0.05, fontface="italic")+
#     xlim_tree(100)
#   
# })
# 
# tree.plot<-plot_grid(plotlist = tree.plots,nrow = 1,labels = names(fitMIXED))
# ggsave("tree_fits.pdf",tree.plot,width = 49,height = 20)
