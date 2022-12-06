setwd("~/PrimateMolars_6traits/EXPLORING_PCM/")
library(ape)
library(phytools)
library(PCMFit)
library(PCMBase)
library(data.table)
library(magrittr)
library(foreach)
library(doParallel)
library(plyr)
library(treeplyr)
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



#Mixed model----
load("~/PrimateMolars_6traits/MixedModels/BIC_results/fitmixed6.RData")
# names(fitMIXED)<-paste("Mixed",1:length(fitMIXED))
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


bestFitMixed<-RetrieveBestFitScore(fitMIXED$`fitMixed 5`)
phy<-bestFitMixed$tree
phy$tip.label<-tree$tip.label
phy<-groupClade(phy,.node = names(phy$part.regime))

regTab<-TipLabelTable(bestFitMixed$inferredModel, offset = 1)
regTab$part.model[]<-"BM"

tree1.plot <- 
  ggtree(phy,size=0,layout = "fan")+
  geom_tree(aes(color=group))+
  scale_color_manual(values=brewer.pal(3,name = "Dark2")[-2])+
  theme(legend.position = "none")
  


bestFitMixed<-RetrieveBestFitScore(fitMIXED$`fitMixed 4`)
phy<-bestFitMixed$tree
phy$tip.label<-tree$tip.label
phy<-groupClade(phy,.node = names(phy$part.regime))

regTab<-TipLabelTable(bestFitMixed$inferredModel, offset = 1)
regTab$part.model[]<-"BM"

tree2.plot <- 
  ggtree(phy,size=0,layout = "fan")+
  geom_tree(aes(color=group))+
  scale_color_manual(values=brewer.pal(3,name = "Dark2"), 
                     labels = c("Ancestral", "Strepsirrhini", "Simiiformes"),
                     name="Regimes")

plot_grid(tree1.plot,tree2.plot,rel_widths = c(1,1.4)) %>% ggsave(., filename = "trees_dists.pdf",width = 12,height = 5)

# tr.plot<-
#   ggplot(SimTr %>% mutate(.,.id=factor(.id, 1:10)), aes(.id, V1))+
#   geom_boxplot(aes(fill=.id), show.legend = F)+
#   # scale_y_log10()+
#   xlab("model")+
#   ylab("Trace")+
#   theme_minimal()+
#   scale_y_log10()

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
  llply(fitMIXED[1:sum(!laply(fitMIXED,is.null))], function(x){
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
ggsave("tree_fits.pdf",tree.plot,width = 49,height = 20)


#PCA plots------
library(mclust)
library(ggforce)

x<-means %>% subset(., M3>0 & M3<1.8) %>% select(.,m1bl:m3md) %>% na.omit %>% log
pca<-prcomp(x)
d_clust <- Mclust(x, G=1:10, 
                  modelNames = mclust.options("emModelNames"))

bics<-data.frame(BIC=d_clust$BIC[,d_clust$modelName],variables="linear", model=d_clust$modelName)

df.pca<-data.frame(x,pca$x,group=factor(d_clust$classification))

# df.means<-data.frame(scale(t(d_clust$parameters$mean),scale=F) %*% pca$rotation, group=factor(1:dim(d_clust$parameters$variance$sigma)[3]))

distPCA.plot<-
  ggplot(df.pca, aes(PC1, PC2))+
  geom_point(aes(color=group), alpha=0.5)+
  # geom_point(aes(color=group), df.means, size=3)+
  stat_ellipse(aes(fill=group),type = "norm",level = 0.5, alpha=0.3,geom = "polygon")+
  # scale_color_manual(values=brewer.pal(3,name = "Dark2")[c(1,3,2)])+
  xlab(paste0("PC1 (",round(100*pca$sdev[1]^2/sum(pca$sdev^2), 2),"%)"))+
  ylab(paste0("PC2 (",round(100*pca$sdev[2]^2/sum(pca$sdev^2), 2),"%)"))+
  theme(legend.position = "none")


x<-means %>% subset(., M3>0 & M3<1.8) %>% select(.,m1a:m3a) %>% na.omit %>% log
pca<-prcomp(x)
d_clust <- Mclust(x, G=1:10, 
                  modelNames = mclust.options("emModelNames"))

bics<-rbind(bics,data.frame(BIC=d_clust$BIC[,d_clust$modelName],variables="area", model=d_clust$modelName))

df.pca<-data.frame(x,pca$x,group=factor(d_clust$classification))

areaPCA.plot<-
  ggplot(df.pca, aes(PC1, PC2))+
  geom_point(aes(color=group))+
  stat_ellipse(aes(fill=group),type = "norm",level = 0.5, alpha=0.3,geom = "polygon")+
  xlab(paste0("PC1 (",round(100*pca$sdev[1]^2/sum(pca$sdev^2), 2),"%)"))+
  ylab(paste0("PC2 (",round(100*pca$sdev[2]^2/sum(pca$sdev^2), 2),"%)"))+
  theme(legend.position = "none")


x<-select(means,M2,M3) %>% na.omit %>% subset(., M3>0 & M3<1.8)
pca<-prcomp(x)
d_clust <- Mclust(x, G=1:10, 
                  modelNames = mclust.options("emModelNames"))
bics<-rbind(bics,data.frame(BIC=d_clust$BIC[,d_clust$modelName],variables="ICM", model=d_clust$modelName))
df.pca<-data.frame(x,pca$x,group=factor(d_clust$classification))

icmPCA.plot<-
  ggplot(df.pca, aes(PC1, PC2))+
  geom_point(aes(color=group))+
  stat_ellipse(aes(fill=group),type = "norm",level = 0.5, alpha=0.3,geom = "polygon")+
  xlab(paste0("PC1 (",round(100*pca$sdev[1]^2/sum(pca$sdev^2), 2),"%)"))+
  ylab(paste0("PC2 (",round(100*pca$sdev[2]^2/sum(pca$sdev^2), 2),"%)"))+
  theme(legend.position = "none")

plot_grid(distPCA.plot+theme(legend.position = "none"),
          areaPCA.plot+theme(legend.position = "none"),
          icmPCA.plot+theme(legend.position = "none"),
          ncol = 1, align = "h",labels = "AUTO") %>%
  ggsave("PCAs.pdf",.,width = 7,height = 7)

mclustBICs.plot<-ggplot(data.frame(i=factor(1:10),bics) %>% mutate(., variables=factor(variables, levels = unique(variables))), aes(i,-BIC))+
  facet_grid(variables~., scale="free")+
  geom_point(aes(shape=model), size=4)+
  geom_line(aes(group=variables))+
  xlab("Cluster number")+
  ylab("BIC")
ggsave("mclustBICs.pdf",mclustBICs.plot,width = 7,height = 7)
  


#----------------









