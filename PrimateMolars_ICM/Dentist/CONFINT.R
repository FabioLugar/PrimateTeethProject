#HEADER------
setwd("~/PrimateMolars/Dentist/")
library(PCMBase)
library(PCMkappa)
library(PCMFit)
library(dentist)
library(ggplot2)
library(cowplot)
library(foreach)
library(plyr)
library(abind)
library(dplyr)
library(expm)
library(xtable)
library(psych)
library(reshape)
options(PCMBase.Threshold.EV = 1e-8)

load("../GlobalModels/fitGlobal.Rdata")

#BM-------

bestModel<-RetrieveBestModel(fitBM)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0","sigma1","sigma2", "sigma3")
wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

best_neglnL<- wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.01,0.01,0.01)
maxPars[]<-c(1.2,1.6,0.5 ,0.5 ,0.5 )
dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

xtable(t(dented_results$all_ranges),digits = 3)

#BMkappa-------

bestModel<-RetrieveBestModel(fitBMkappaG)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0","kappa")
best_neglnL<- wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.1)
maxPars[]<-c(1.2,1.6,0.9)
dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

xtable(t(dented_results$all_ranges),digits = 3)

#OU_kappaG-------
bestModel<-RetrieveBestModel(fitOUkappaG)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0","kappa", "H1", "H2", "H3","theta1","theta2")
best_neglnL<- c(-fitOUkappaG$logLikOptim)
# PCMLik(X,tree, bestModel, SEsp, metaI=PCMBaseCpp::PCMInfoCpp)

wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.1,0.01,-0.5,0.01,0.8,0.8)
maxPars[]<-c(1.2,1.6,0.9,0.5 , 0.5,0.5 ,1.2,1.6)

dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)


xtable(t(dented_results$all_ranges),digits = 3)

pairwiseParameter.plots<-
  combn(names(best_par),2) %>%
  alply(., 2, function(x){
    dent_results2plot<-
      dented_results$results[dented_results$acceptances,] %>% 
      mutate(., small=(neglnL<(dented_results$best_neglnL+2))) %>%
      subset(., neglnL<dented_results$best_neglnL+10) %>%
      select(., any_of(c(x[1],x[2],"small"))) %>%
      arrange(., small)
    colnames(dent_results2plot)[-3]<-c("trait1","trait2")
    
    dent_results2plot
    
    ggplot(dent_results2plot, aes(trait1, trait2))+
      geom_point(aes(fill=small), shape=21, show.legend = T, show.legend=F)+
      stat_ellipse(aes(trait1, trait2),subset(dent_results2plot, small==T))+
      scale_fill_manual(values=c("white", "black"))+
      xlab(x[1])+
      ylab(x[2])
  })

pdf("pairwiseParameter.pdf",width = 7,height = 7)
pairwiseParameter.plots
dev.off()


confintML<-dented_results$results %>%
  subset(., neglnL<(dented_results$best_neglnL+2)) %>% unique

# mod<-bestModel
# 
# stat_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
#   PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
#   stationary <-PCMkappa:::StationaryVariance(mod$H[,,1],mod$Sigma_x[,,1])
#   statElli<-ellipse::ellipse(stationary, centre=mod$X0)
#   data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
# }
# 
# sigma_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
#   PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
#   statElli<-ellipse::ellipse(tcrossprod(mod$Sigma_x[,,1]), centre=mod$Theta[,1])
#   data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
# }
# 
# 
# selec_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
#   PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
#   H<-PCMApplyTransformation(mod)$H[,,1]
#   Sigma<-tcrossprod(mod$Sigma_x[,,1])
#   statElli<-ellipse::ellipse(
#     solve(sqrtm(H))%*%Sigma%*%solve(sqrtm(H)),centre=mod$Theta[,1])
#   data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
# }
# names(selec_elip)[-1]<-c("x","y","x.1","y.1")
# 
# meansvars<-
#   means %>% select(., Species,M2, M3) %>% na.omit()
# meansvars<-
#   cbind(meansvars,
#         M2sd=M2errors[meansvars$Species,"sd"],
#         M3sd=M3errors[meansvars$Species,"sd"]
#   ) %>% 
#   mutate(., 
#          M2_max= M2+1.96*M2sd,
#          M2_min= M2-1.96*M2sd,
#          M3_max= M3+1.96*M3sd,
#          M3_min= M3-1.96*M3sd)
# 
# 
# stat_ellip.plot<- 
#   rbind(data.frame(sigma_elip,matrix="'A.Genetic Drift ' (Myr^{-1})"),
#         data.frame(selec_elip,matrix="'B.Adaptive landscape'")) %>%#,
#         # data.frame(stat_elip,matrix="'C.Stationary variance'")) %>%
#   ggplot(.)+
#   facet_grid(.~matrix, labeller = label_parsed)+
#   geom_pointrange(aes(x=M2,y=M3,ymin=M3_min, ymax=M3_max), 
#                   meansvars,fill="white",shape=21)+
#   geom_pointrange(aes(x=M2,y=M3,xmin=M2_min, xmax=M2_max), 
#                   meansvars,fill="white",shape=21)+
#   # geom_point(aes(x=x,y=y, group=.id, color=matrix), alpha=0.1, show.legend = F)+
#   geom_segment(aes(x=x,y=y,xend=x.1,yend=y.1, group=.id, color=matrix), alpha=0.1, show.legend = F)+
#   theme_minimal()+
#   coord_cartesian()+
#   xlab("M2/M1")+
#   ylab("M3/M1")+
#   coord_fixed()+
#   ylim(c(0,2.3))+
#   scale_color_manual(values=c("deeppink", "darkviolet","darkseagreen4" ))
# stat_ellip.plot
# ggsave("stat_ellip2.pdf",stat_ellip.plot,width = 18, height = 11, units="cm")


H<-fitOU$modelOptim$H

hxs<-confintML %>%
  ddply(., 1, function(x) {
    h<-matrix(0,2,2)
    s<-tcrossprod(bestModel$Sigma_x[,,1])*x[,"kappa"]
    H[,,1][upper.tri(s,diag = T)]<-unlist(x[,c("H1","H2","H3")])
    H<-PCMApplyTransformation(H)
    data.frame(sigma=tr(s),H=tr(H[,,1]), deltaL=abs(best_neglnL-x[,"neglnL"]))
  })
# library(lmodel2)

hxs.plot<-
  ggplot(hxs, aes(sigma, H))+
  geom_point()+
  geom_smooth(method="lm",se=T, color="black")+
  xlab(expression(tr(Sigma)))+
  ylab(expression(tr(H)))

ggsave("hxs.pdf",hxs.plot,width = 5, height = 5)

#OU-------

bestModel<-RetrieveBestModel(fitOU)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0", "H1", "H2", "H3","theta1","theta2","sigma1","sigma2", "sigma3")
best_neglnL<- c(-fitOU$logLikOptim)
# PCMLik(X,tree, bestModel, SEsp, metaI=PCMBaseCpp::PCMInfoCpp)

wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.01,-0.5,0.01,0.8,0.8,0.01,0.01,0.01)
maxPars[]<-c(1.2,1.6,0.5 , 0.5,0.5 ,1.2,1.6,0.5 ,0.5 ,0.5 )

dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

xtable(t(dented_results$all_ranges),digits = 3)


dented_results$all_ranges[,c("sigma1","sigma2", "sigma3")]

H<-fitOU$modelOptim$H

sigmaG<-PCMApplyTransformation(fitOUkappaG$modelOptim)$Sigma_x[,,1] %>% tcrossprod

hxs<-dented_results$results %>% subset(., neglnL<(best_neglnL+8)) %>% unique %>%
  ddply(., 1, function(x) {
    s<-h<-matrix(0,2,2)
    s[upper.tri(s,diag = T)]<-unlist(x[,c("sigma1","sigma2","sigma3")])
    s<-tcrossprod(s)
    H[,,1][upper.tri(s,diag = T)]<-unlist(x[,c("H1","H2","H3")])
    H<-PCMApplyTransformation(H)
    data.frame(sigma=tr(s),H=tr(H[,,1]), deltaL=abs(best_neglnL-x[,"neglnL"]))
  })

#OUdiag--------

bestModel<-RetrieveBestModel(fitOUdiag)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0", "H1", "H2","theta1","theta2","sigma1","sigma2", "sigma3")
best_neglnL<- c(-fitOUdiag$logLikOptim)
# PCMLik(X,tree, bestModel, SEsp, metaI=PCMBaseCpp::PCMInfoCpp)

wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.01,0.01,0.8,0.8,0.01,0.01,0.01)
maxPars[]<-c(1.2,1.6,0.5 ,0.5 ,1.2,1.6,0.5 ,0.5 ,0.5 )

dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

xtable(t(dented_results$all_ranges),digits = 3)


dented_results$all_ranges[,c("sigma1","sigma2", "sigma3")]

H<-fitOU$modelOptim$H

sigmaG<-PCMApplyTransformation(fitOUkappaG$modelOptim)$Sigma_x[,,1] %>% tcrossprod

hxs<-dented_results$results %>% subset(., neglnL<(best_neglnL+8)) %>% unique %>%
  ddply(., 1, function(x) {
    s<-h<-matrix(0,2,2)
    s[upper.tri(s,diag = T)]<-unlist(x[,c("sigma1","sigma2","sigma3")])
    s<-tcrossprod(s)
    H[,,1][upper.tri(s,diag = T)]<-unlist(x[,c("H1","H2","H3")])
    H<-PCMApplyTransformation(H)
    data.frame(sigma=tr(s),H=tr(H[,,1]), deltaL=abs(best_neglnL-x[,"neglnL"]))
  })

#OUkappaGdiag--------

bestModel<-RetrieveBestModel(fitOUkappaGdiag)
best_par<-PCMParamGetShortVector(bestModel)
names(best_par)<-c("X0","Y0","kappa", "H1", "H2","theta1","theta2")
best_neglnL<- c(-fitOUkappaGdiag$logLikOptim)
# PCMLik(X,tree, bestModel, SEsp, metaI=PCMBaseCpp::PCMInfoCpp)

wPCMLik<-function(par, model, SE=SEsp, tree=tree, X=X){
  PCMParamLoadOrStore(model, par, offset = 0 , load=T)
  -c(PCMLik(X,tree,model,SE, metaI=PCMBaseCpp::PCMInfoCpp))
}

wPCMLik(best_par,bestModel,SEsp,tree,X)

minPars<-maxPars<-best_par
minPars[]<-c(0.8,0.8,0.1,0.01,0.01,0.8,0.8)
maxPars[]<-c(1.2,1.6,0.9,0.5 ,0.5 ,1.2,1.6)

set.seed(999)
dented_results <- dent_walk(par=best_par, 
                            fn=wPCMLik, 
                            best_neglnL=best_neglnL,
                            nsteps=10000, 
                            print_freq=250,lower_bound = minPars,upper_bound = maxPars,
                            X=X, tree=tree, model=bestModel,SE = SEsp)

xtable(t(dented_results$all_ranges),digits = 3)

confintML<-dented_results$results %>%
  subset(., neglnL<(dented_results$best_neglnL+2)) %>% unique

mod<-bestModel

stat_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
  PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
  stationary <-PCMkappa:::StationaryVariance(mod$H[,,1],mod$Sigma_x[,,1])
  statElli<-ellipse::ellipse(stationary, centre=mod$X0)
  data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
}

sigma_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
  PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
  statElli<-ellipse::ellipse(tcrossprod(mod$Sigma_x[,,1]), centre=mod$Theta[,1])
  data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
}


selec_elip<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
  PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
  H<-PCMApplyTransformation(mod)$H[,,1]
  Sigma<-tcrossprod(mod$Sigma_x[,,1])
  statElli<-ellipse::ellipse(
    (solve(sqrtm(H))%*%Ps$G%*%solve(sqrtm(H)))-Ps$pooled,
    centre=mod$Theta[,1])
  data.frame(.id=i,statElli[-nrow(statElli),],statElli[-1,])
}
names(selec_elip)[-1]<-c("x","y","x.1","y.1")

meansvars<-
  means %>% select(., Species,M2, M3) %>% na.omit()
meansvars<-
  cbind(meansvars,
        M2sd=M2errors[meansvars$Species,"sd"],
        M3sd=M3errors[meansvars$Species,"sd"]
  ) %>% 
  mutate(., 
         M2_max= M2+1.96*M2sd,
         M2_min= M2-1.96*M2sd,
         M3_max= M3+1.96*M3sd,
         M3_min= M3-1.96*M3sd)


stat_ellip.plot<- 
  rbind(data.frame(sigma_elip,matrix="'A.Genetic Drift ' (Myr^{-1})"),
        data.frame(selec_elip,matrix="'B.Adaptive landscape'")) %>%#,
  # data.frame(stat_elip,matrix="'C.Stationary variance'")) %>%
  ggplot(.)+
  facet_grid(.~matrix, labeller = label_parsed)+
  geom_pointrange(aes(x=M2,y=M3,ymin=M3_min, ymax=M3_max), 
                  meansvars,fill="white",shape=21)+
  geom_pointrange(aes(x=M2,y=M3,xmin=M2_min, xmax=M2_max), 
                  meansvars,fill="white",shape=21)+
  # geom_point(aes(x=x,y=y, group=.id, color=matrix), alpha=0.1, show.legend = F)+
  geom_segment(aes(x=x,y=y,xend=x.1,yend=y.1, group=.id, color=matrix), alpha=0.05, show.legend = F)+
  theme_minimal()+
  coord_cartesian()+
  xlab("m2/m1")+
  ylab("m3/m1")+
  coord_fixed()+
  ylim(c(0,2.3))+
  scale_color_manual(values=c("deeppink", "darkviolet","darkseagreen4" ))
stat_ellip.plot
ggsave("stat_ellip2.pdf",stat_ellip.plot,width = 18, height = 11, units="cm")


ICM_regr<-foreach(i=1:dim(confintML)[1],.combine = "rbind") %do% {
  PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
  stationary <-PCMkappa:::StationaryVariance(mod$H[,,1],mod$Sigma_x[,,1])
  optima<-mod$Theta[,1]
  slope=eigen(stationary)$vector[2,1]*1/eigen(stationary)$vector[1,1]
  intercept=optima[2]-slope*optima[1]
  data.frame(slope,intercept)
  }


regr.plot<-
  ggplot()+
  geom_abline(aes(slope=slope, intercept=intercept), ICM_regr,alpha=0.1)+
  # geom_point(aes(x=M2,y=M3), meansvars,fill="white",shape=21)+
  geom_pointrange(aes(x=M2,y=M3,ymin=M3_min, ymax=M3_max),
                  meansvars,fill="white",shape=21)+
  geom_pointrange(aes(x=M2,y=M3,xmin=M2_min, xmax=M2_max),
                  meansvars,fill="white",shape=21)+
  geom_abline(slope=2, intercept=-1, linetype=2, color="magenta")+
  theme_minimal()+
  coord_cartesian()+
  xlab("m2/m1")+
  ylab("m3/m1")+
  # coord_fixed()+
  ylim(c(0,2.3))

ggsave("regr.pdf",regr.plot,width = 7, height = 5)

t1.2<-
  foreach(i=1:dim(confintML)[1],.combine = "rbind") %dopar% {
    PCMParamLoadOrStore(mod, unlist(confintML[i,-1]), offset = 0 , load=T)
    # stationary <-PCMkappa:::StationaryVariance(mod$H[,,1],mod$Sigma_x[,,1])
    modt<-PCMApplyTransformation(mod)
    H<-modt$H[,,1]
    Sigma<-tcrossprod(modt$Sigma_x[,,1])
    Omega<-solve(sqrtm(H))%*%Sigma%*%solve(sqrtm(H))

    optima<-mod$Theta[,1]
    
    V<-eigen(Omega)$vector
    Hr<-t(V)%*%H%*%V
    data.frame(M2=diag(H)[1],
               M3=diag(H)[2],
               A.I=diag(Hr)[1],
               dev=diag(Hr)[2]) %>%
      mutate_all(.,.funs = function(x) log(2)/x)
  }

t12.plot<-
  t1.2 %>% melt %>% 
  ggplot(., aes(variable, value))+
  geom_violin(aes(fill=variable), show.legend = F)+
  xlab("")+
  ylab(expression(paste("Phylogenetic Half-life (",t[1/2],")",sep = "")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_x_discrete(
                   labels=c("m2/m1","m3/m1", "Activation-Inhibition gradient", "Deviation from ICM"))

ggsave("t12.pdf", t12.plot, width = 7, height = 5)
