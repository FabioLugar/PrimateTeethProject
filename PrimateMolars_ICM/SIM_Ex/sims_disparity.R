setwd("~/PrimateMolars/SIM_Ex/")
library(ape)
library(phytools)
library(PCMFit)
library(PCMBaseCpp)
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
doParallel::registerDoParallel(cores = 40)
mcaffinity(20:60)
load("~/PrimateMolars/GlobalModels/fitGlobal.Rdata")

mod<-RetrieveBestModel(fitOUkappaG) %>% PCMApplyTransformation
simsOU<-foreach(i=seq_len(1000)) %dopar% {
  PCMSim(tree=tree_na, model=mod, X0=mod$X0, SE=SEsp[,tree_na$tip.label])
}

mod$H[,,]<-0
simsBM<-foreach(i=seq_len(1000)) %dopar% {
  PCMSim(tree=tree_na, model=mod, X0=mod$X0, SE=SEsp[,tree_na$tip.label])
}

ellis<-rbind(
  foreach(i=seq_along(simsOU),.combine = "rbind") %do% {
  x<-simsOU[[i]]
  statElli<-ellipse::ellipse(cov(t(x)), centre=mod$X0)
  data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,],type="OU")
  },
  foreach(i=seq_along(simsBM),.combine = "rbind") %do% {
  x<-simsBM[[i]]
  statElli<-ellipse::ellipse(cov(t(x)), centre=mod$X0)
  data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,],type="BM")
})



disparity_sim.plot<-
  ggplot(ellis %>% subset(.,.id %in% seq_len(100)), aes(x,y))+
  geom_segment(aes(x=x,y=y,xend=x.1,yend=y.1, group=.id, color=type), alpha=1, show.legend = F)+
  facet_grid(.~type)+
  coord_fixed()+
  geom_point(aes(M2,M3),data.frame(t(X)), shape=21,fill="white")+
  theme_bw()+
  coord_fixed()+
  xlab("m2/m1")+
  ylab("m3/m1")+
  ylim(c(0.25,2.3))
ggsave("disparity_sim.pdf",disparity_sim.plot,width = 8,height = 4)


disparity.plot<-
  data.frame(OU=simsOU %>% laply(., function(x) tr(cov(t(x[,1:Ntip(tree_na)])))),
           BM=simsBM %>% laply(., function(x) tr(cov(t(x[,1:Ntip(tree_na)]))))) %>%
  melt %>%
  ggplot(., aes(variable, sqrt(value)))+
  geom_violin(aes(fill=variable), show.legend = FALSE, alpha=0.5)+
  ylab("Disparity")+
  xlab("")+
  scale_y_log10()+
  scale_fill_manual(values=c("darkviolet","darkseagreen4" ))+
  geom_hline(yintercept = sqrt(tr(cov(X_na))))


physig.sims<-
  data.frame(OU=simsOU %>% 
               laply(., function(x) physignal(t(x[,1:Ntip(tree_na)]), tree_na,iter = 0)$phy.signal, .parallel = T),
           BM=simsBM %>% 
             laply(., function(x) physignal(t(x[,1:Ntip(tree_na)]), tree_na,iter = 0)$phy.signal, .parallel = T))

phylosig.plot<-
  physig.sims[,2:1] %>%
  melt %>%
  ggplot(., aes(variable, value))+
  theme_bw()+
  geom_violin(aes(fill=variable), show.legend = FALSE, alpha=0.5)+
  ylab("Phylogenetic Signal (K)")+
  xlab("")+
  scale_y_log10()

phylosig.plot


ggsave("phylosig.pdf",phylosig.plot,width = 5,height = 5)
ggsave("disparity.pdf",disparity.plot,width = 5,height = 5)

both.plot<-plot_grid(disparity_sim.plot, phylosig.plot,labels = "AUTO",ncol = 1,align = "v")

ggsave("disparity_phylosig.pdf",both.plot,width = 10,height = 10)

save.image("sims.Rdata")

tr(var(na.omit(t(X))))
tr(var(t(X), use="complete.obs"))

disparity.plot+geom_hline(yintercept = tr(var(t(X), use="complete.obs")))


