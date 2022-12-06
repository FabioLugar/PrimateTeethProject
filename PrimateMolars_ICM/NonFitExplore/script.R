setwd("~/PrimateMolars/NonFitExplore/")
library(mvtnorm)
library(emorph2)
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
library(treeplyr)
library(phytools)
library(reshape2)
library(scales)
library(PCMkappa)
library(PCMFit)
library(phyloEM)
library(evolqg)
library(dispRity)
library(paleotree)

load("../MixedModels/BIC_results/fitMIXED.RData")

bestModel<-RetrieveBestModel(fitOUkappaG)

stationary <-PCMkappa:::StationaryVariance(bestModel$H[,,1],bestModel$Sigma_x[,,])

PCAphylo<-eigen(stationary)

pcascore<-  scale(na.omit(t(X)), scale=F) %*% PCAphylo$vectors

ages<-setNames(node.depth.edgelength(tree_na),c(tree_na$tip.label,tree_na$node.label))
pcascore<-data.frame(actin=c(pcascore[,1],fastAnc(tree_na,pcascore[,1])),
                     dev=c(pcascore[,2],fastAnc(tree_na,pcascore[,2])),
                     ages)

opt<-data.frame(opt=t((bestModel$Theta[,]-colMeans(na.omit(t(X)))) %*% PCAphylo$vectors),
                sd=sqrt(PCAphylo$values),
                PC=c("Activation-inhibition","Deviation from ICM"))

TT.plot<-
  data.frame(rbindlist(list(
  data.frame(pcascore[tree_na$edge[,1],c(1,3)],end=pcascore[tree_na$edge[,2],c(1,3)], PC="Activation-inhibition"),
  data.frame(pcascore[tree_na$edge[,1],c(2,3)],end=pcascore[tree_na$edge[,2],c(2,3)], PC="Deviation from ICM")))) %>%
  ggplot(.)+
  facet_grid(PC~.,scales = "free")+
  # geom_smooth(aes(ages,actin),method = "loess",linetype=0, fill="red", alpha=0.7)+
  geom_hline(aes(yintercept=opt), opt, size=1, color="goldenrod",alpha=1)+
  geom_point(aes(x=75, y=opt), opt, size=5, color="goldenrod",alpha=1)+
  geom_segment(aes(x=ages,y=actin,xend=end.ages,yend=end.actin), alpha=0.5)+
  geom_point(aes(end.ages, end.actin), size=0.5)+
  geom_point(aes(ages, actin), size=0.5)+
  theme_bw()+
  xlab("age")+
  ylab("")
  
ggsave("TT.pdf", TT.plot,width = 8,height = 7)

ltt_obj<-ltt(tree)

wd<-5

interp<-
  foreach(i=1:nrow(tree_na$edge),.combine = "rbind") %do% {
    nods<-tree_na$edge[i,]
    spls<-floor(pcascore[nods[2],"ages"]-pcascore[nods[1],"ages"])
    if(spls>0){
      data.frame(pc1=seq(pcascore[nods[1],1],pcascore[nods[2],1],length.out = spls+2),
                 pc2=seq(pcascore[nods[1],2],pcascore[nods[2],2],length.out = spls+2),
                 ages=seq(pcascore[nods[1],3],pcascore[nods[2],3],length.out = spls+2))
    }
}
rvar<-foreach(i=0:round(max(interp$ages)-wd),.combine = "rbind") %do%{
  x<-interp
  data.frame(ages=i+(wd/2), var=psych::tr(var(subset(x, ages>i & ages<(i+wd))[,1:2])))
}

var.plot<-
  ggplot(rbind(data.frame(rvar, type="disparity"),
               data.frame(ages=ltt_obj$times, var=log(ltt_obj$ltt)*0.02, type="lineages")), 
         aes(ages-max(ages), var))+
  geom_line(aes(linetype=type))+
  # geom_point()+
  theme_bw()+
  xlab("age")+
  ylab("Variance")+
  scale_y_continuous(
    name = "Disparity",
    sec.axis = sec_axis(~./0.02, name="log-lineages"))+
  scale_linetype_manual(name="", values=c(2,1))
var.plot
ggsave("var.pdf", var.plot,width = 8,height = 3.5)

simsvars<-foreach(i=seq_len(1000),.combine = "cbind") %dopar%{
  Xs<-PCMSim(tree=tree_na, model=bestModel, X0=bestModel$X0, SE=SEsp[,tree_na$tip.label])
  pcascore<-  scale(na.omit(t(Xs)), scale=F) %*% PCAphylo$vectors
  pcascore<-data.frame(actin=pcascore[,1],
                       dev=pcascore[,2],
                       ages)
  
  interp<-
    foreach(i=1:nrow(tree_na$edge),.combine = "rbind") %do% {
      nods<-tree_na$edge[i,]
      spls<-floor(pcascore[nods[2],"ages"]-pcascore[nods[1],"ages"])
      if(spls>0){
        data.frame(pc1=seq(pcascore[nods[1],1],pcascore[nods[2],1],length.out = spls+2),
                   pc2=seq(pcascore[nods[1],2],pcascore[nods[2],2],length.out = spls+2),
                   ages=seq(pcascore[nods[1],3],pcascore[nods[2],3],length.out = spls+2))
      }
    }
  foreach(i=0:round(max(interp$ages)-wd),.combine = "c") %do%{
    x<-interp
    psych::tr(var(subset(x, ages>i & ages<(i+wd))[,1:2]))
  }
}



var_confint.plot<-
  ggplot(data.frame(rvar, conf=t(apply(simsvars,1,quantile, prob=c(0.025,0.975)))), 
         aes(ages-max(ages), var))+
  geom_line(aes())+
  # geom_point()+
  theme_bw()+
  xlab("age")+
  ylab("Disparity")+
  geom_ribbon(aes(ymin=conf.2.5., ymax=conf.97.5.), alpha=0.3)
var_confint.plot
ggsave("var_confint.pdf", var_confint.plot,width = 8,height = 3.5)


library(data.table)
library(MASS)

a<-data.table(rvar)

a[,merge:=ages]

b<-data.table(ages=ltt_obj$times, lineages=log(ltt_obj$ltt))

b[,merge:=ages]

setkeyv(a,c('merge'))

setkeyv(b,c('merge'))

Merged=b[a,roll='nearest']

roblm<-rlm(var ~ lineages, data = Merged)

ggplot(Merged, aes(lineages, var))+
  geom_abline(intercept=roblm$coef[1], slope=roblm$coef[2])+
  geom_point()+
  theme_bw()

roblm$residuals 

ggplot(cbind(Merged[,.(ages,lineages=lineages*0.022,var)]) %>% melt(., id.vars=c("ages")), 
       aes(ages-max(ages), value))+
  geom_line(aes(linetype=variable))+
  geom_point(aes(ages-max(ages), var,fill=res), cbind(Merged[,.(ages,var)],res=roblm$residuals), shape=21, size=2)+
  theme_bw()+
  xlab("age")+
  ylab("Variance")+
  scale_y_continuous(
    name = "Disparity",
    sec.axis = sec_axis(~./0.022, name="log-lineages"))+
  scale_linetype_manual(name="", values=c(2,1))+
  scale_fill_gradient(low="black",high = "white")

revts(ggtree(tree) + theme_tree2())

# data.frame(evol=Evolvability(Ps$G),flex=Flexibility(Ps$G)) %>%
#   melt() %>%
#   ggplot(., aes(variable, value))+
#   geom_violin()+
#   facet_wrap(.~variable,scales = "free")+
#   geom_point(aes(variable, value), data.frame(evol=Evolvability(Ps$G,beta.mat = matrix(c(1,0,0,1, PCAphylo$vectors),2)),flex=Flexibility(Ps$G,beta.mat = matrix(c(1,0,0,1, PCAphylo$vectors),2)),
#                           name=c("M2","M3","PC1", "PC2")
#   ) %>% melt)+
#   geom_label_repel(aes(variable, value, label=name), data.frame(evol=Evolvability(Ps$G,beta.mat = matrix(c(1,0,0,1, PCAphylo$vectors),2)),flex=Flexibility(Ps$G,beta.mat = matrix(c(1,0,0,1, PCAphylo$vectors),2)),
#                           name=c("M2","M3","PC1", "PC2")
#   ) %>% melt)
# 
# b<-apply(rmvnorm(1000,sigma = diag(2)),1,Normalize)

Gdist<-llply(1:1000,function(x){
  x<-rmvnorm(100, mean = Gmeans,sigma = Gcov)
  x<-x[,1:3]*x[,4:6]
  x<-x[,2:3]/x[,1]
  BayesianCalculateMatrix(lm(x~1),samples = 2)$Ps[1,,]
}) %>% abind(., along=3) %>% aperm(., c(3,1,2))


b<-matrix(c(1,0,0,1, PCAphylo$vectors),2)
colnames(b)<-c("M2","M3","PC1", "PC2")

evol<-adply(Gdist,1,Evolvability,beta.mat = b)[,-1] 
evol_minmax<-t(adply(Gdist,1,Evolvability)[,-1] %>% apply(., 1, range))
evol<-melt((evol-evol_minmax[,1])/evol_minmax[,2])

sp <- split(evol$value, evol$variable)
a <- lapply(seq_along(sp), function(i){
  d <- density(sp[[i]],adjust = 2)
  k <- which.max(d$y)
  data.frame(value = names(sp)[i], xmax = d$x[k], ymax = d$y[k])
})
a <- do.call(rbind, a)

evol.plot<-
  ggplot(evol, aes(value))+
  geom_density(aes(fill=variable), alpha=0.5, show.legend = F,adjust = 2)+
  theme_bw()+
  geom_label(data = a, 
            aes(x = xmax, y = ymax, 
                label = value, vjust = -0.5))+
  ylim(c(0,45))+
  scale_fill_brewer(palette = "Set1")+
  ylab("")+
  xlab("Evolvability")

ggsave("evol.pdf", evol.plot,width = 5,height = 5)

