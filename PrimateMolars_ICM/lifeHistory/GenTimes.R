setwd("~/PrimateMolars/lifeHistory/")
library(readxl)
library(plyr)
library(dplyr)
library(ape)
library(geiger)
library(phytools)
# AnAge<-read.csv("anage_data.csv") %>% subset(., Order=="Primates") %>% mutate(., GSP=paste(Genus,Species, sep="_"))
AnAge<-read.csv("annasdata/anagedata.csv") %>% subset(., Order=="Primates") %>% mutate(., GSP=paste(Genus,Species, sep="_"))
AnAge$GSP
tree<-read.nexus("../../WD/UpdatedData/medianTree.nexus") %>% ladderize()
means<-read.csv("../../WD/UpdatedData/means.csv")
# means<-read.csv("annasdata/means.csv")

subNames<-read.csv("../../data/trees/subNamesmeta.csv") 
#apply name changes to original dataset
for(i in 1:dim(subNames)[1]){
  AnAge$GSP<-sub(subNames[i,1], subNames[i,2], AnAge$GSP)
}

means<-subset(means, Species %in% tree$tip.label)
sum(AnAge$GSP %in% means$Species)
AnAge<-AnAge[AnAge$GSP %in% tree$tip.label,]

gentimes<-setNames(AnAge$Female.maturity..days.+AnAge$Gestation.Incubation..days.,AnAge$GSP) %>% na.omit
gentimes<-gentimes/365.25

treeData<-treedata(tree,data.frame(gentimes))
treeData$phy
treeData$data

anc<-fastAnc(treeData$phy, treeData$data[,1])

1000000/anc[1]
