setwd("~/PrimateMolars_6traits/Bootstrap/")
library(ape)
library(phytools)
library(PCMFit)
library(PCMBase)
library(PCMkappa)
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

load("~/PrimateMolars_6traits/MixedModels/BIC_results/fitmixed6.RData")

numBootstraps <- 100
bestFitMixed<-RetrieveBestFitScore(fitMIXED[[4]])
bestModel<-bestFitMixed$inferredModel
startree<-stree(length(tree$tip.label), "star")
startree$edge.length<-rep(max(nodeHeights(tree)), times=length(tree$tip.label))
X0<-bestModel$X0


SimTr<-ldply(bestModel[-1], function(x){
  model<-PCM(sub("Omitted_X0","Global_X0",class(x)[1]), k = 6)
  model$X0<-X0
  for(i in names(x)) model[i]<-x[i]
  ldply(1:100, function(i) {
    XX<-PCMSim(startree, model, X0=X0)
    tr(var(t(XX)))
  },.parallel = T)
})

tr.plot<-
  ggplot(SimTr, aes(.id, V1))+
  geom_boxplot(aes(fill=.id), show.legend = F)+
  # scale_y_log10()+
  xlab("model")+
  ylab("Trace")+
  theme_minimal()

tr.plot
saveRDS(SimTr,"SimTr.RDS")


SimTr6t<-readRDS("~/PrimateMolars_6traits/Bootstrap/SimTr.RDS")
SimTr3t<-readRDS("~/PrimateMolars_3traits/Bootstrap/SimTr.RDS")
SimTr2t<-readRDS("~/PrimateMolars/Bootstrap/SimTr.RDS")
SimTr3t$.id<-factor(SimTr3t$.id)
levels(SimTr3t$.id)<-letters[1:3]
SimTr6t$.id<-factor(SimTr6t$.id)
levels(SimTr6t$.id)<-letters[1:3]

sims.plot<-
  rbind(data.frame(traits="ICM",SimTr2t),
      data.frame(traits="Area",SimTr3t),
      data.frame(traits="Distances",SimTr6t)) %>% 
  mutate(., .id=revalue(.id,  c(a = "Ancestral", b= "Strepsirrhini", c = "Simiiformes"))) %>%
  ggplot(., aes(.id,V1))+
  facet_grid(traits~.,scales = "free")+
  geom_violin(aes(fill=.id), show.legend = F, alpha=0.7)+
  xlab("")+
  ylab("Disparity")+
  scale_fill_brewer(type="qual", palette = 2)

ggsave("sims.pdf",sims.plot,width = 6,height = 6)
