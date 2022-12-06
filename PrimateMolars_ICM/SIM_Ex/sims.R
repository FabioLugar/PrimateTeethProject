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
doParallel::registerDoParallel(cores = 50)

load("~/PrimateMolars/GlobalModels/fitGlobal.Rdata")

mod<-RetrieveBestModel(fitOUkappaG)
sims<-foreach(i=seq_len(1000)) %dopar% {
  PCMSim(tree=tree, model=mod, X0=mod$X0, SE=SEsp)
}


ages<-setNames(node.depth.edgelength(tree),c(tree$tip.label,tree$node.label))
ages<-ages[1:(Ntip(tree))] #getting only tips
ages<-ages[ages<max(ages)-0.5] #removing extant data

longtips<-sort(unique(tree$edge[tree$edge.length>1,2])) #getting tips with leading long edges
longtips<-longtips[longtips<Ntip(tree)] #geting only tips
longtips<-tree$tip.label[longtips] #tips with long leading edges
ex_sp<-longtips[longtips %in% names(ages)] #getting only extinct branches
ex_sp<-ex_sp[!ex_sp %in% tree$tip.label[is.na(colSums(X))]]


obsD<-mahalanobis(t(X[,ex_sp]), c(mod$Theta),Ps$G) #empirical mahalanobis distance from 
simD<-ldply(sims, function(X) mahalanobis(t(X[,ex_sp]), c(mod$Theta),Ps$G)) #simulated mahalanobis distance from 

rbind(data.frame(type="sim",D=unlist(simD)),
      data.frame(type="obs",D=unlist(obsD))) %>%
  ggplot(., aes(type, D^2))+
  # geom_boxplot(fill="gray")+
  geom_violin(fill="gray")+
  # geom_jitter()+
  scale_y_log10()+
  theme_bw()


D<-aaply(as.matrix(simD),1,median)^2 %>% log
x.dens <- density(D)
df.dens <- data.frame(x = x.dens$x, y = x.dens$y)
perc<-quantile(D, c(0.025, 0.975))
densityD.plot<-
  data.frame(D) %>%
  ggplot(., aes(D))+
  theme_bw()+
  geom_area(data = subset(df.dens, x >= perc[1] & x <= perc[2]), 
            aes(x=x,y=y), fill = 'gray')+
  geom_density(aes(x=D, y = ..density..), color = 'black')+
  geom_vline(xintercept=median(obsD)^2 %>% log)+
  xlab(expression(D^2))+
  ylab("Density")
densityD.plot
ggsave("densityD.pdf",densityD.plot,width = 5,height = 5)

confInt<-foreach(i=ex_sp,.combine = "rbind") %do% {
  data.frame(sp=i, ages=ages[i], sup=quantile(simD[,i],probs = c(0.95)))
}


  
densityDtime.plot<-
  data.frame(confInt,D=obsD) %>%
  ggplot(., aes(ages, D))+
  theme_bw()+
  geom_ribbon(aes(x= ages, y=NULL,ymin = 1, ymax = sup), alpha=0.3)+
  geom_point()
densityDtime.plot
ggsave("densityDtime.pdf",densityDtime.plot,width = 7,height = 5)

