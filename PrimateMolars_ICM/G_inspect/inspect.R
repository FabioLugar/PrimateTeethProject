setwd("~/PrimateMolars/G_inspect/")
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
library(ggplot2)
library(RColorBrewer)
library(emorph2)
library(mvtnorm)
doParallel::registerDoParallel(cores = 50)
options(PCMBase.Threshold.EV = 1e-8)
load("../GlobalModels/fitGlobal.Rdata")

covA<-diag(c(emorph2:::varProd(Gcov[c(1,4),c(1,4)],Gmeans[c(1,4)]),
             emorph2:::varProd(Gcov[c(2,5),c(2,5)],Gmeans[c(2,5)]),
             emorph2:::varProd(Gcov[c(3,6),c(3,6)],Gmeans[c(3,6)])))
covA[1,2] <- covA[2,1] <-emorph2:::covProd(Gcov[c(1,2,4,5),c(1,2,4,5)],Gmeans[c(1,2,4,5)])
covA[1,3] <- covA[3,1] <-emorph2:::covProd(Gcov[c(1,3,4,6),c(1,3,4,6)],Gmeans[c(1,3,4,6)])
covA[2,3] <- covA[3,2] <-emorph2:::covProd(Gcov[c(2,3,5,6),c(2,3,5,6)],Gmeans[c(2,3,5,6)])
meansA<-Gmeans[1:3]*Gmeans[4:6]

Ecovs<-foreach(i=1:10000) %dopar%{
  E<-rmvnorm(300,sigma =  Gcov, mean = Gmeans)
  cov(E[,c(1:3)]*E[,c(4:6)])
}

cors<-laply(Ecovs, function(M) cor(M[upper.tri(M, diag=T)], covA[upper.tri(covA, diag=T)]),.parallel = T)
difs<-ldply(Ecovs, function(M) {
  covA<-covA / meansA%*%t(meansA)
  M   <-M / meansA%*%t(meansA)
  M[upper.tri(M, diag=T)]-covA[upper.tri(covA, diag=T)]
  },.parallel = T)

ggplot(data.frame(cor=cors),aes(cor))+ 
  geom_density(fill="black",alpha=0.4)+
  xlim(c(0,1))

# names(difs)<-c("sigma[M1]^2","sigma[M1,M2]","sigma[M2]^2","sigma[M1,M3]","sigma[M2,M3]","sigma[M3]^2")
difs<-difs[,c(1,3,6,2,4,5)]

names(difs)<-c("sigma[M1]^2","sigma[M2]^2","sigma[M3]^2","sigma[M1M2]","sigma[M1M3]","sigma[M2M3]")
error.plot<-ggplot(melt(difs), aes(variable, value))+
                geom_violin(draw_quantiles = c(0.025,0.975))+
                scale_x_discrete(labels=parse(text=c(names(difs))))+
                xlab("")+
                ylab("Difference (Coefficient of Variation)")+
                ylim(c(-0.003,0.003))+
                geom_hline(yintercept = 0, linetype=2)
          
ggsave("error.pdf",error.plot,width = 6,height = 4)
