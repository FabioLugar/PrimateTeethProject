#HEADER-----
setwd("~/PrimateMolars/Rates/")
library(ggtree)
# library(RRphylo)
library(emorph2)
library(data.table)
library(ape)
# library(phylobase)
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
library(PCMFit)
library(mvtnorm)
library(evolqg)
library(scales)

load("../sp_means_trees_updated.Rdata")
means_red<-select(means,Species,M2,M3) %>% 
  # subset(., M3>0) %>% 
  na.omit
# class(tree)<-class(tree)[-1]
matcheddat <- make.treedata(tree, means_red)
matcheddat$dat<-as.data.frame(matcheddat$dat)
rownames(matcheddat$dat)<-matcheddat$phy$tip.label

# estimating gentimes----

AnAge<-read.csv("../lifeHistory/anage_data.csv") %>% subset(., Order=="Primates") %>% mutate(., GSP=paste(Genus,Species, sep="_"))

subNames<-read.csv("../subNamesmeta.csv") 
#apply name changes to original dataset
for(i in 1:dim(subNames)[1]){
  AnAge$GSP<-sub(subNames[i,1], subNames[i,2], AnAge$GSP)
}
AnAge<-AnAge[AnAge$GSP %in% tree$tip.label,]
gentimes<-setNames(AnAge$Female.maturity..days.+AnAge$Gestation.Incubation..days.,AnAge$GSP) %>% na.omit
gentimes<-gentimes/365.25
treeData<-make.treedata(tree,data.frame(gentimes))
anc<-fastAnc(treeData$phy, setNames(unlist(treeData$dat),treeData$phy$tip.label),CI = TRUE)
gt<-anc$CI95[1,]

#getting scaling factor-----
ndz<-length(tree$edge.length)
tc<-1000000
tg<-round(tc/ceiling(gt))

tgs<-runif(n = 1000,tg[2],tg[1])
Nes<-runif(n = 1000,20000,1000000)

tab4imcov<-dplyr::select(molarMeasures %>% subset(., M3!=0), Species, M2:M3) %>% na.omit
tab4imcov<-data.table(tab4imcov)[, if(.N > 30) .SD,Species]

#distributions for G, expected D and LGGD-----
# Gdist<-BayesianCalculateMatrix(lm(cbind(M2,M3)~Species,tab4imcov),samples = 1000)$Ps
# for(i in 1:1000) Gdist[i,,]<- Gdist[i,,] * runif(1,0.3,0.6)
load("../Gbaboon2.Rdata")
Gdist<-llply(1:1000,function(x){
  x<-rmvnorm(100, mean = Gmeans,sigma = Gcov)
  x<-x[,1:3]*x[,4:6]
  x<-x[,2:3]/x[,1]
  BayesianCalculateMatrix(lm(x~1),samples = 2)$Ps[1,,]
}) %>% abind(., along=3) %>% aperm(., c(3,1,2))

NULLDexp<-aaply(1:100, 1, function(i) {
  Gdist[i,,] * tgs[i]/Nes[i]
})
NULLlandeDist<-aaply(1:100, 1, function(i) {
  Dexp<-Gdist[i,,] * tgs[i]/Nes[i]
  Dzs<-rmvnorm(ndz,sigma = Dexp)
  mahalanobis(Dzs, center = FALSE, cov = Gdist[i,,]) * Nes[i]/tgs[i]
})
CI<-quantile(melt(NULLlandeDist)$value,probs = c(0.025,0.975))
matcheddat$phy<-multi2di(matcheddat$phy)
dzs<-apply(matcheddat$dat,2,pic,phy=matcheddat$phy)
lengths<-pic(matcheddat$dat$M2,phy=matcheddat$phy,var.contrasts = T)[,2]

#calculate "empirical" LGGD----
landeDists<-aaply(1:100, 1, function(i) {
  mahalanobis(dzs, center = FALSE, cov = Gdist[i,,]) * Nes[i]/tgs[i]
})

confint95.stat <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.975))
  names(r) <- c("ymin","ymax")
  r
}

rate_age.plot<-
  landeDists %>% t %>% data.frame(.,age=aaply(as.numeric(colnames(landeDists)), 1, nodeheight, tree=matcheddat$phy)) %>%
  melt(., id.vars=c("age")) %>%
  ggplot(., aes(age-max(age),value))+
  stat_summary(fun.data=confint95.stat, geom = "errorbar", width=.1, alpha=0.5)+
  stat_summary(fun="median",size = 2, geom = "point")+
  # geom_point(aes(age-max(age), m, fill=length),data.frame(m=apply(landeDists, 2, median),
                              # length=lengths,
                              # age=aaply(as.numeric(colnames(landeDists)), 1, nodeheight, tree=matcheddat$phy)),pch=21)+
  # geom_jitter()+
  scale_y_log10()+
  scale_x_continuous(breaks=c(-60,-40,-20,0),labels=c(60,40,20,0))+
  geom_hline(aes(yintercept=CI),data.frame(CI=CI), linetype=2)+
  xlab("Age (myr)")+
  ylab("LGGD")+
  theme_bw()#+
  # scale_fill_gradient(low="#f4f1de", high="#cc2936",trans="log",breaks=c(1,7,54))
test<-cor.test(apply(landeDists, 2, median),lengths, method="spearman")

rate_length.plot<-ggplot(data.frame(m=apply(landeDists, 2, median),length=lengths), aes(length,m))+
  geom_smooth(method="lm")+geom_point()+
  scale_y_log10()+scale_x_log10()+
  ylab("LGGD")+
  xlab("Branch length (Myr)")+
  annotate("text",x=10,y=9e+01,
           label=as.expression(bquote("Spearman's "~rho~"="~.(round(test$estimate,3))~", p<1e-18")))



# ------------

phy<-matcheddat$phy
# g1 = as(phy, 'phylo4')
# d = data.frame(color=sample(NA, phy$Nnode+1, replace=T))
# rownames(d) = phy$tip.label
# g2 = phylo4d(g1)

rateLong<-rep(NA, length(phy$tip.label)+phy$Nnode)
for(i in unique(phy$edge[,1])){
  j<-phy$edge[,1]==i
  rateLong[phy$edge[j,2]]<-colMeans(log(landeDists))[i-length(phy$tip.label)]
}
rateLong[is.na(rateLong)]<-colMeans(log(landeDists))[1]

# rNodeData <- data.frame(rate = exp(rateLong),
#                         row.names = 1:length(rateLong))
# tdata(g2) <- rNodeData



rNodeData <- data.frame(node = 1:(length(phy$tip.label)+phy$Nnode),
                        trait = exp(rateLong),
                        row.names = c(phy$tip.label,
                                      length(phy$tip.label)+(1:phy$Nnode)))


phy<-full_join(phy, rNodeData, by = 'node')
# 
# 
#   geom_tree(,phy, continuous = 'colour', size=1, show.legend=F)

tree_LGGD.plot<-
  ggtree(phy, aes(color=trait),size=2,layout="fan")+
  scale_color_gradientn(colors=c("darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","gray","gray","gray","gray","mediumpurple4","mediumpurple4","mediumpurple4","mediumpurple4"),
  trans="log", breaks=round(setNames(CI,NULL),2), name="LGGD")+
  theme(legend.position = c(0.8, 0.65))+
  geom_tiplab(aes(label=sub("_"," ",label)), size=1, hjust=-0.05, fontface="italic")+
  xlim_tree(100)

ggsave("rate_age.pdf",rate_age.plot,width = 11,height = 10, units="cm")
ggsave("rate_length.pdf",rate_length.plot,width = 14,height = 10, units="cm")
ggsave("tree_LGGD.pdf",tree_LGGD.plot,width = 10,height = 10)

