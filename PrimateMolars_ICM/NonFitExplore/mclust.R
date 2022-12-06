setwd("~/PrimateMolars/NonFitExplore/")
load("../sp_means_trees_updated.Rdata")
library(mclust)
library(dplyr)

d_clust.ICM <- select(means,M2,M3) %>% na.omit %>% subset(., M3>0 & M3<1.8) %>%
  Mclust(., G=1:10, 
         modelNames = mclust.options("emModelNames"))

d_clust.ICM$BIC
plot(d_clust.ICM)
4
0
pca<-means %>% subset(., M3>0 & M3<1.8) %>% select(.,m1a:m3a) %>% na.omit %>% log %>% prcomp
d_clust.area <- 
  Mclust(pca$x, G=1:10, 
         modelNames = mclust.options("emModelNames"))

d_clust.area$BIC
plot(d_clust.area)
2
0

pca<-means %>% subset(., M3>0 & M3<1.8) %>% select(.,m1bl:m3md) %>% na.omit %>% log %>% prcomp

d_clust.dist <- 
  Mclust(pca$x[,], G=1:10, 
         modelNames = mclust.options("emModelNames"))

d_clust.dist$BIC
plot(d_clust.dist)
