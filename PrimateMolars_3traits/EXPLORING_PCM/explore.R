setwd("~/PrimateMolars_3traits/EXPLORING_PCM/")
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
# library(ggtree)
library(RColorBrewer)
library(psych)
library(geiger)
library(mvtnorm)
library(ggstance)
library(evolqg)
library(reshape)
library(expm)
options(PCMBase.Threshold.EV = 1e-8)

# doParallel::registerDoParallel(cores = 50)
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


#Mixed model----
load("~/PrimateMolars_3traits/MixedModels/BIC_results/fitMIXED.RData")
# load("~/PrimateMolars_3traits/MixedModels/BIC_results/fitmixed6.RData")

names(fitMIXED)<-paste("Mixed",1:length(fitMIXED))

names(fitMIXED)<-paste("fitMixed",1:10)
modelCompare<-
  rbind(list(BM=fitBM,
             OU=fitOU,
             BMkappaP=fitBMkappaP,
             OUkappaP=fitOUkappaP,
             BMkappaG=fitBMkappaG,
             OUkappaG=fitOUkappaG) %>%
          ldply(.,function(x) data.frame(logLik=x$logLikOptim,
                                         n=dim(x$X)[2],
                                         k=PCMParamCount(x$modelOptim))),
        ldply(fitMIXED[!laply(fitMIXED,is.null)],function(x){
          mod<-RetrieveBestFitScore(x)$inferredModel
          data.frame(logLik=x$tableFitsRR$logLik[2],
                     n=x$tableFitsRR$nobs[2],
                     k=PCMParamCount(mod))
        })) %>%
  mutate(
    .,
    AICc = (2 * k - 2 * logLik) + (2 * k) * (k + 1) / (n - k - 1),
    BIC = (log(n) * k - 2 * logLik)
  )
arrange(modelCompare, BIC)

library(xtable)
xtable(modelCompare[,c(1,4,2,5)])


write.csv(modelCompare,"modelCompare.csv")

bestFitMixed<-RetrieveBestFitScore(fitMIXED$`fitMixed 2`)
# 
# startree<-stree(length(tree$tip.label), "star")
# startree$edge.length<-rep(max(nodeHeights(tree)), times=length(tree$tip.label))
# X0<-bestFitMixed$inferredModel$X0
# 
# SimTr<-ldply(bestFitMixed$inferredModel[-1], function(x){
#   model<-PCM(sub("Omitted_X0","Global_X0",class(x)[1]), k = 3)
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

# attr(bestFitMixed$inferredModel, "tree")$tip.label<-tree$tip.label
# 
# bestFitMixed<-RetrieveBestFitScore(fitMIXED$`fitMixed 2`)
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
#   theme_minimal()+
#   scale_y_log10()
# 
# trait.plot<-
#   PCMPlotTraitData2D(attr(bestFitMixed$inferredModel, "X")[c(1,3),],
#                    attr(bestFitMixed$inferredModel, "tree"),
#                    alphaPoints = 1,scaleSizeWithTime = F,
#                    numTimeFacets = 1)
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
# trait.plot<-ggplot(df2plot, aes(m1a, m2a))+
#   geom_point(aes(x=m1a,y=m2a, color=regime))+
#   # geom_smooth(aes(color=regime),method="lm",se=F)+
#   theme_bw()
# 
# ggsave("Mixed_model_vars.pdf",plot_grid(tree.plot,trait.plot,tr.plot,ncol = 2),width = 10,height = 10)


#########


tree.plots<-
  llply(fitMIXED[!laply(fitMIXED,is.null)], function(x){
  bestFit <- RetrieveBestFitScore(x)
  modelTree<-bestFit$tree
  
  modelTree$tip.label<-tree$tip.label
  PCMTreePlot(modelTree)%<+% 
    TipLabelTable(bestFit$inferredModel, offset = 5) + 
    # geom_tiplab2(size = 2) + 
    geom_tiplab(aes(label = part.model), offset = 25) + 
    geom_tiplab(aes(label=sub("_"," ",label)), size=2, hjust=-0.05, fontface="italic")+
    xlim_tree(100)
  
})

tree.plot<-plot_grid(plotlist = tree.plots,nrow = 1,labels = names(fitMIXED))
tree.plot
ggsave("tree_fits.pdf",tree.plot,width = 49,height = 20)

# subset(means, Genus=="Hylobates") %>% select(., m1a,m2a,m3a) %>% mutate(., m1a=log(m1a),m2a=log(m2a),m3a=log(m3a)) %>% apply(.,2,sd)
# 
# cld<-c("Sphacorhysis_burntforkensis", "Walshina_shifrae", "Walshina_esmaraldensis", "Walshina_mcgrewi", "Trogolemur_fragilis", "Anthradapis_vietnamensis", "Nosmips_aenigmaticus", "Nannopithex_filholi", "Nannopithex_zuccolae", "Melaneremia_bryanti", "Tetonoides_pearcei", "Anemorhysis_pattersoni", "Arapahovius_gazini", "Melaneremia_schrevei", "Arapahovius_advena", "Anemorhysis_natronensis", "Anemorhysis_sublettensis", "Anemorhysis_wortmani", "Teilhardina_tenuicula", "Anemorhysis_savagei", "Teilhardina_crassidens")
# subset(means, Species %in% cld) %>% select(., m1a,m2a,m3a) %>% mutate(., m1a=log(m1a),m2a=log(m2a),m3a=log(m3a)) %>% apply(.,2,sd, na.rm=T)