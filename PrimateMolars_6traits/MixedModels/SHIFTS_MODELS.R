setwd("~/PrimateMolars_6traits/MixedModels/")
library(ggtree)
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
load("../GlobalModels/fitGlobal.Rdata")
# source('DefineParameterLimits.R', local=FALSE)
doParallel::registerDoParallel(cores = 30)
# openblasctl::openblas_set_num_threads(5)
parallel::mcaffinity(seq(1:30))
options(PCMBase.Threshold.EV = 1e-8)
rm(i)

generatePCMModelsFunction <- function() {
  PCMGenerateModelTypes()
  fileName <- '../GlobalModels/DefineParameterLimits_GLOBAL.R'
  codeDefineLimits <- readChar(fileName, file.info(fileName)$size)
  eval(parse(text = codeDefineLimits), .GlobalEnv)
}

tax2fix<-head(names(which(table(means$tax)>5)),-1)
clades_fixed<-alply(tax2fix,1,function(x) {
  sps<-subset(subset(means, Species %in% tree$tip.label), tax==x)$Species
  node<-getMRCA(tree,sps)
  extract.clade(tree, node)$node.label
}) %>% unlist %>% sort
other_clades<-c("791", "508","651", "797", "795", "790","652", "647")

# (481, 794, 493, 595)
pdf(file = "skiped_clades.pdf",height = 10,width = 7)
ggtree(tree) + 
  geom_point2(aes(subset=(node %in% clades_fixed)),color="blue",size=2)+
  geom_point2(aes(subset=(node %in% other_clades)),color="red",size=2)
dev.off()

clades_fixed<-c(clades_fixed,other_clades)

nruns<-10
fitMIXED<-vector("list",nruns)
for(j in 1:length(fitMIXED)) {
  
  #Fitting unknown shifts-----
  
  prefixFiles = paste0("mixed_",j)
  
  currentResultFile <- paste0("Current_", prefixFiles, ".RData")
  if(file.exists(currentResultFile)) {
    load(currentResultFile)
    tableFitsPrev <- listResults$tableFits
  } else {
    tableFitsPrev <- NULL
  }
  
  
  fitMIXED[[j]] <- PCMFitMixed(modelTypes = MGPMDefaultModelTypes()[c("B","E")],
                               subModels = NULL,
                               X = X, tree = tree, metaIFun = PCMBaseCpp::PCMInfoCpp,
                               generatePCMModelsFun = generatePCMModelsFunction,
                               maxNumRoundRobins = 1, maxNumPartitionsInRoundRobins = 2,
                               tableFitsPrev = tableFitsPrev,
                               prefixFiles = prefixFiles,
                               minCladeSizes = 35L,
                               skipNodes = clades_fixed,
                               SE=SEsp,
                               doParallel = TRUE,
                               scoreFun = BIC)
  
  save.image(paste0("BIC_results/fitmixed", j, ".RData"))
}
unlink(paste0("BIC_results/fitmixed", 1:nruns, ".RData"))
save.image("BIC_results/fitMIXED.RData")

