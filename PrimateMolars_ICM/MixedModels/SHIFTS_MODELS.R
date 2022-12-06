setwd("~/PrimateMolars/MixedModels/")
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
doParallel::registerDoParallel(cores = 50)
options(PCMBase.Threshold.EV = 1e-8)
load("../GlobalModels/fitGlobal.Rdata")
rm(i)

generatePCMModelsFunction <- function() {
  PCMGenerateModelTypes()
  fileName <- 'DefineParameterLimits_shifts.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

nruns<-10
fitMIXED<-vector("list",nruns)
for(j in 1:length(fitMIXED)) {
  
  #Fitting unknown shifts-----
  
  prefixFiles = paste0("mixed_",j)
  
  currentResultFile <- paste0("Current_",j,"_", prefixFiles, ".RData")
  if(file.exists(currentResultFile)) {
    load(currentResultFile)
    tableFitsPrev <- listResults$tableFits
  } else {
    tableFitsPrev <- NULL
  }
  
  # tax2fix<-head(names(which(table(means$tax)>5)),-1)
  # clades_fixed<-alply(tax2fix,1,function(x) {
  #   sps<-subset(subset(means, Species %in% tree$tip.label), tax==x)$Species
  #   node<-getMRCA(tree,sps)
  #   extract.clade(tree, node)$node.label
  # }) %>% unlist %>% sort
  
  fitMIXED[[j]] <- PCMFitMixed(modelTypes = MGPMDefaultModelTypes()[c("B","E")],
                               subModels = NULL,
                               X = X, tree = tree, metaIFun = PCMBaseCpp::PCMInfoCpp,
                               generatePCMModelsFun = generatePCMModelsFunction,
                               maxNumRoundRobins = 1, maxNumPartitionsInRoundRobins = 2,
                               tableFitsPrev = tableFitsPrev,
                               prefixFiles = prefixFiles,
                               minCladeSizes = 5L,
                               # skipNodes = clades_fixed,
                               SE=SEsp,
                               doParallel = TRUE,
                               scoreFun = BIC)
  
  save.image(paste0("fitmixed", j, ".RData"))
}
unlink(paste0("fitmixed", 1:nruns, ".RData"))
save.image("fitMIXED.RData")
