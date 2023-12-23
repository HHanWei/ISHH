################################################################################
# This document includes codes for conducting simulation studies.
################################################################################
rm(list=ls())
gc()
library(Matrix)
library(MASS)
library(glmnet)
source('main-function.R')
source('sim-function.R')
################################################################################
n <- 120; q <- 3; p <- 12; s <- 4                     
rho0 <- 0.3; rho1 <- 0.3; rho2 <- 0.7                    
pr1 <- c(0.25,0.25,0.25,0.25); pr2 <- c(0.2,0.2,0.3,0.3)      
mu1 <- 1; mu2 <- 2; va <- 0.25                    
sp1 <- c(1,2,3,4)             
sp2 <- c(1,2,3,5)             
sp3 <- c(1,2,5,6)             
sp4 <- c(1,5,6,7)  
################################################################################
# Example
# mu1-AR1 with balanced design under prior1
wholedata <- Generatedata(n, q, p, rho0, rho1, pr1, s, mu1, va, banded = F)
init <- initFun(wholedata, lambda.sparse = seq(0.01,1,length=50), div.num = 4, lambda1 = 1, lambda2 = 0.5)
para <- list(lambda1=seq(0.21,0.01,length=5), lambda2=seq(0.85,0.65,length=5))
tau <- c(1,seq(0.9,0,length=10))
tunelambda <- tuninglambda(wholedata, init$xi.init, init$gamma.init, para, epsi = 0.01)
prior1 <- ADMM0(wholedata, init$xi.init, init$gamma.init,
                lambda1 = (tunelambda$per.df.tune)$lambda1, lambda2 = (tunelambda$per.df.tune)$lambda2, 
                sp = sp1, admmstep = 1, iter.max = 50, epsi = 0.01, mcp.para = 3, penal.para = 1)
tunesp1 <- tuningtau(wholedata, init$xi.init, init$gamma.init, xi.prior = prior1$xi.gf, gamma.prior = prior1$gamma.gf, 
                     sp = sp1, tau = tau, lambda1.opt = (tunelambda$per.df.tune)$lambda1, lambda2.opt = (tunelambda$per.df.tune)$lambda2)




