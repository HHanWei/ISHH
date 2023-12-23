################################################################################
# This document includes codes for conducting real data analysis.
################################################################################
############################### Part 0 #########################################
rm(list=ls())
gc()
library(Matrix)
library(MASS)
library(glmnet)
source('main-function.R')
source('case-function.R')
################################################################################
load("LUADdata116.RData")
# data processing ##############################################################
# Find rough group structure ###################################################
freq <- data$freq
gene.ind <- which(freq>=50)
sp1 <- which(freq[gene.ind,]>=300)
sp2 <- which(freq[gene.ind,]>=400)
sp3 <- which(freq[gene.ind,]>=500)
sp4 <- which(freq[gene.ind,]>=800)
lambda2 <- seq(0.01,0.4,length=10)
res <- vector(mode = 'list', length = length(lambda2))
bic <- c()
for (i in 1:length(lambda2)) {
  res[[i]] <- casepreFun(x.mat = as.matrix(data$data.x.scale), y.mat = as.matrix(data$data.y.scale), lambda1 = 0.001, lambda2 = lambda2[i], epsi = 0.001)
  bic[i] <- res[[i]]$BIC.var 
}
res.final <- res[[min(which(bic==min(bic)))]]
group.rough <- res.final$group.rough
wholedata <- list(data.x = as.matrix(data$data.x.scale), data.z = as.matrix(data$data.z.scale[,gene.ind]), data.y = as.matrix(data$data.y.scale),
                  group.rough = group.rough, sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4)
save(wholedata, file='wholedata.RData')
################################################################################

################################################################################
rm(list=ls())
gc()
library(Matrix)
library(MASS)
library(glmnet)
################################################################################
source('main-function.R')
source('case-function.R')
source('sim-function.R')

############################### Part I #########################################
load('wholedata.RData')
# Our proposed method
init <- caseinitFun(wholedata, lambda.sparse = seq(0.01,0.05,length = 10), div.num = c(2,3), lambda1 = 0.001, lambda2 = 0.001)
para <- list(lambda1=seq(0.30,0.10,length=9), lambda2=seq(0.45,0.15,length=13))
tau <- c(1,seq(0.9,0,length=10))
tunelambda <- tuninglambda(wholedata, init$xi.init, init$gamma.init, para, epsi = 0.001, multi = T)
# scenario 1
prior1 <- ADMM0(wholedata, init$xi.init, init$gamma.init,
                lambda1 = (tunelambda$per.df.tune)$lambda1, lambda2 = (tunelambda$per.df.tune)$lambda2, 
                sp = wholedata$sp1, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
tunesp1 <- tuningtau(wholedata, init$xi.init, init$gamma.init, prior1$xi.gf, prior1$gamma.gf, 
                     sp = wholedata$sp1, tau = tau, lambda1.opt = (tunelambda$per.df.tune)$lambda1, lambda2.opt = (tunelambda$per.df.tune)$lambda2, epsi = 0.001, multi = T)
# scenario 2
prior2 <- ADMM0(wholedata, init$xi.init, init$gamma.init,
                lambda1 = (tunelambda$per.df.tune)$lambda1, lambda2 = (tunelambda$per.df.tune)$lambda2, 
                sp = wholedata$sp2, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
tunesp2 <- tuningtau(wholedata, init$xi.init, init$gamma.init, prior2$xi.gf, prior2$gamma.gf, 
                     sp = wholedata$sp2, tau = tau, lambda1.opt = (tunelambda$per.df.tune)$lambda1, lambda2.opt = (tunelambda$per.df.tune)$lambda2, epsi = 0.001, multi = T)
# scenario 3
prior3 <- ADMM0(wholedata, init$xi.init, init$gamma.init,
                lambda1 = (tunelambda$per.df.tune)$lambda1, lambda2 = (tunelambda$per.df.tune)$lambda2, 
                sp = wholedata$sp3, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
tunesp3 <- tuningtau(wholedata, init$xi.init, init$gamma.init, prior3$xi.gf, prior3$gamma.gf, 
                     sp = wholedata$sp3, tau = tau, lambda1.opt = (tunelambda$per.df.tune)$lambda1, lambda2.opt = (tunelambda$per.df.tune)$lambda2, epsi = 0.001, multi = T)
# scenario 4
prior4 <- ADMM0(wholedata, init$xi.init, init$gamma.init,
                lambda1 = (tunelambda$per.df.tune)$lambda1, lambda2 = (tunelambda$per.df.tune)$lambda2, 
                sp = wholedata$sp4, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
tunesp4 <- tuningtau(wholedata, init$xi.init, init$gamma.init, prior4$xi.gf, prior4$gamma.gf, 
                     sp = wholedata$sp4, tau = tau, lambda1.opt = (tunelambda$per.df.tune)$lambda1, lambda2.opt = (tunelambda$per.df.tune)$lambda2, epsi = 0.001, multi = T)
# summary
tunesp <- list(tunesp1 = tunesp1, tunesp2 = tunesp2, tunesp3 = tunesp3, tunesp4 = tunesp4)
################################################################################

############################### Part II ########################################
# Use ANOVA to check the significance of subgroups on some clinical features
# load("wholedata.RData")
# load("tunesp.RData")
group.fine <- tunesp$tunesp1$res.tune$gr.gf
id <- rownames(wholedata$data.x)
group.label <- rep(0,116)
for (i in 1:4) {
  group.label[group.fine[[i]]] <- rep(i,length(group.fine[[i]]))
}

all <- read.csv('LUAD_clinical.csv')
all$PATIENT_ID <- gsub('\\-','\\.',all$PATIENT_ID)
all116 <- all[match(id,all$PATIENT_ID),]
rownames(all116) <- NULL
all116 <- cbind(data.frame('group.label'=group.label),all116)
for (i in 1:nrow(all116)) {
  for (j in 1:ncol(all116)) {
    if(all116[i,j] == '[Not Available]'){
      all116[i,j] <- NA
    }
  }
}

aov.data1 <- na.omit(all116[,c('group.label','PATIENT_ID','CARBON_MONOXIDE_DIFFUSION_DLCO')])
x <- aov(CARBON_MONOXIDE_DIFFUSION_DLCO~group.label, data = aov.data1)
summary(x) #0.0586

aov.data2 <- na.omit(all116[,c('group.label','PATIENT_ID','FEV1_FVC_RATIO_POSTBRONCHOLIATOR')])
x <- aov(FEV1_FVC_RATIO_POSTBRONCHOLIATOR~group.label, data = aov.data2)
summary(x) #0.0139

aov.data3 <- na.omit(all116[,c('group.label','PATIENT_ID','FEV1_FVC_RATIO_PREBRONCHOLIATOR')])
x <- aov(FEV1_FVC_RATIO_PREBRONCHOLIATOR~group.label, data = aov.data3)
summary(x) #0.0666
################################################################################

################################ Part III ######################################
# Remove 1 data point
# load("wholedata.RData")
index <- function(x, y){
  index <- c()
  for (i in 1:length(x)) {
    index[i] <- which(y==x[i])
  }
  return(index)
}
# Left 115 data points
wholedata115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  wholedata115[[i]] <- wholedata
  wholedata115[[i]]$data.x <- wholedata$data.x[-i,]
  wholedata115[[i]]$data.z <- wholedata$data.z[-i,]
  wholedata115[[i]]$data.y <- as.matrix(wholedata$data.y[-i,])
  wholedata115[[i]]$group.rough <- list(index(setdiff(wholedata$group.rough[[1]],i),setdiff(c(1:116),i)),
                                        index(setdiff(wholedata$group.rough[[2]],i),setdiff(c(1:116),i)))
}
# Our proposed method
init115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  init115[[i]] <- caseinitFun(wholedata115[[i]], lambda.sparse = seq(0.01,0.05,length = 10), div.num = c(2,3), lambda1 = 0.001, lambda2 = 0.001)
}
# scenario 1
prior1115 <- vector(mode = 'list', length = 116)
tunesp1115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  cat('-----------', i, '-th data' ,'--------------\n')
  prior1115[[i]] <- ADMM0(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, lambda1 = 0.275, lambda2 = 0.45, 
                          sp = wholedata115[[i]]$sp1, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
  tunesp1115[[i]] <- tuningtau(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, 
                               xi.prior = prior1115[[i]]$xi.gf, gamma.prior = prior1115[[i]]$gamma.gf, sp = wholedata115[[i]]$sp1, 
                               tau = 0.6, lambda1.opt = 0.275, lambda2.opt = 0.45, epsi = 0.001, multi = T)
}
# scenario 2
prior2115 <- vector(mode = 'list', length = 116)
tunesp2115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  cat('-----------', i, '-th data' ,'--------------\n')
  prior2115[[i]] <- ADMM0(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, lambda1 = 0.275, lambda2 = 0.45, 
                          sp = wholedata115[[i]]$sp2, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
  tunesp2115[[i]] <- tuningtau(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, 
                               xi.prior = prior2115[[i]]$xi.gf, gamma.prior = prior2115[[i]]$gamma.gf, sp = wholedata115[[i]]$sp2, 
                               tau = 0.6, lambda1.opt = 0.275, lambda2.opt = 0.45, epsi = 0.001, multi = T)
}
# scenario 3
prior3115 <- vector(mode = 'list', length = 116)
tunesp3115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  cat('-----------', i, '-th data' ,'--------------\n')
  prior3115[[i]] <- ADMM0(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, lambda1 = 0.275, lambda2 = 0.45, 
                          sp = wholedata115[[i]]$sp3, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
  tunesp3115[[i]] <- tuningtau(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, 
                               xi.prior = prior3115[[i]]$xi.gf, gamma.prior = prior3115[[i]]$gamma.gf, sp = wholedata115[[i]]$sp3, 
                               tau = 0.8, lambda1.opt = 0.275, lambda2.opt = 0.45, epsi = 0.001, multi = T)
}
# scenario 4
prior4115 <- vector(mode = 'list', length = 116)
tunesp4115 <- vector(mode = 'list', length = 116)
for (i in 1:116) {
  cat('-----------', i, '-th data' ,'--------------\n')
  prior4115[[i]] <- ADMM0(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, lambda1 = 0.275, lambda2 = 0.45, 
                          sp = wholedata115[[i]]$sp4, admmstep = 1, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, multi = T)
  tunesp4115[[i]] <- tuningtau(wholedata = wholedata115[[i]], xi.init = init115[[i]]$xi.init, gamma.init = init115[[i]]$gamma.init, 
                               xi.prior = prior4115[[i]]$xi.gf, gamma.prior = prior4115[[i]]$gamma.gf, sp = wholedata115[[i]]$sp4, 
                               tau = 0.8, lambda1.opt = 0.275, lambda2.opt = 0.45, epsi = 0.001, multi = T)
}






