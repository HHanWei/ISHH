#####################################################################################
## All functions for simulation studies.
#####################################################################################

###########################################################################################################################
## Part I #################################################################################################################
##################################### Functions for generating simulation data ############################################
Generatecor <- function(q, p, rho1, rho2, banded = F, corr = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generating the correlation structure
  ## -----------------------------------------------------------------------------------------------------------------
  sig1 <- matrix(rho1, nrow = q, ncol = q)
  diag(sig1) <- 1
  sig2 <- matrix(0, nrow = p, ncol = p)
  if(banded){
    for (i in 1:(p-1)) {
      sig2[i,i+1] <- 0.5
    }
    sig2 <- sig2+t(sig2)
    diag(sig2) <- 1
  }
  else{
    for (i in 1:p) {
      for (j in i:p) {
        sig2[i,j] <- rho2^abs(i-j)
        sig2[j,i] <- sig2[i,j]
      }
    }
  }
  
  if(corr){
    min <- min(min(sig1),min(sig2))
    sig <- matrix(0.05*min, nrow = q+p, ncol = q+p)
    sig[1:q,1:q] <- sig1
    sig[(q+1):(q+p),(q+1):(q+p)] <- sig2
  }
  else{
    sig <- matrix(0, nrow = q+p, ncol = q+p)
    sig[1:q,1:q] <- sig1
    sig[(q+1):(q+p),(q+1):(q+p)] <- sig2
  }
  return(sig)
}

Generategroup <- function(n, pr, Hgroup = list(c(1,2),c(3,4))){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generating the hierarchical subgroup structure 
  ## -----------------------------------------------------------------------------------------------------------------
  group.rough <- vector(mode = 'list', length = length(Hgroup))
  group.true <- vector(mode = 'list', length = length(pr))
  n.num <- sample(n, n, replace = F)
  
  group.true[[1]] <- sort(n.num[1:(n*pr[1])])
  for (k2 in 2:length(pr)) {
    group.true[[k2]] <- sort(n.num[(n*sum(pr[1:(k2-1)])+1):(n*sum(pr[1:k2]))])
  }
  
  for (k1 in 1:length(Hgroup)) {
    Hgroupk1 <- Hgroup[[k1]]
    group.rough[[k1]] <- sort(unlist(group.true[Hgroupk1]))
  }
  
  return(list(group.rough = group.rough, group.true = group.true))
}

Generatecoef <- function(n, q, p, s, mu, group.rough, group.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generating the hierarchical regression coefficients
  ## -----------------------------------------------------------------------------------------------------------------
  xi1 <- matrix(mu, nrow = q, ncol = 1); xi2 <- -xi1
  xi.true.list <- list(xi1, xi2)
  alpha11 <- matrix(0, nrow = p, ncol = 1)
  alpha11[1:s,] <- rep(0.5*mu,s); alpha21 <- -alpha11
  alpha12 <- matrix(0, nrow = p, ncol = 1)
  alpha12[1:s,] <- rep(2*mu,s); alpha22 <- -alpha12
  alpha.true.list <- list(alpha11, alpha12, alpha21, alpha22)
  
  xi.true <- rbind(xi1, xi2)
  alpha.true <- rbind(alpha11, alpha12, alpha21, alpha22)
  
  beta.true <- as.matrix(LMatrixFun(q, group.rough)%*%xi.true)
  gamma.true <- as.matrix(LMatrixFun(p, group.true)%*%alpha.true)
  
  beta.true.mat <- matrix(beta.true, nrow = q, ncol = n)
  gamma.true.mat <- matrix(gamma.true, nrow = p, ncol = n)
  
  beta.true.list <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    beta.true.list[[i]] <- beta.true.mat[,i]
  }
  gamma.true.list <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    gamma.true.list[[i]] <- gamma.true.mat[,i]
  }
  
  return(list(xi.true = xi.true, alpha.true = alpha.true, beta.true = beta.true, gamma.true = gamma.true,
              beta.true.mat = beta.true.mat, gamma.true.mat = gamma.true.mat,
              xi.true.list = xi.true.list, alpha.true.list = alpha.true.list, 
              beta.true.list = beta.true.list, gamma.true.list = gamma.true.list))
}

Generatedata <- function(n, q, p, rho1, rho2, pr, s, mu, va, banded = F, corr = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generating wholedata 
  ## -----------------------------------------------------------------------------------------------------------------
  sig <- Generatecor(q, p, rho1, rho2, banded = banded, corr = corr)
  xz <- MASS::mvrnorm(n, rep(0,(q+p)), sig)
  data.x <- xz[,1:q]
  data.z <- xz[,(q+1):(q+p)]
  
  group <- Generategroup(n, pr)
  group.rough <- group$group.rough
  group.true <- group$group.true
  
  coef <- Generatecoef(n, q, p, s, mu, group.rough, group.true)
  
  beta.true.mat <- coef$beta.true.mat
  gamma.true.mat <- coef$gamma.true.mat
  
  err <- MASS::mvrnorm(n, 0, va)
  
  data.y <- matrix(0, nrow = n, ncol = 1)
  for (i in 1:n) {
    data.y[i,] <- data.x[i,]%*%beta.true.mat[,i]+data.z[i,]%*%gamma.true.mat[,i]+err[i,]
  }
  
  return(list(data.y = data.y, data.x = data.x, data.z = data.z, group.rough = group.rough, group.true = group.true,
              xi.true = coef$xi.true, alpha.true = coef$alpha.true, beta.true = coef$beta.true, gamma.true = coef$gamma.true,
              beta.true.mat = beta.true.mat, gamma.true.mat = gamma.true.mat,
              xi.true.list = coef$xi.true.list, alpha.true.list = coef$alpha.true.list, 
              beta.true.list = coef$beta.true.list, gamma.true.list = coef$gamma.true.list))
}


###########################################################################################################################
## Part II ################################################################################################################
######################################## Functions for evaluating performances  ###########################################
grperFun <- function(group, group.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate estimated group structure performance (RI, TPR, FPR)
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group.true))
  cap.matrix <- pairmatFun(group)
  cap.true.matrix <- pairmatFun(group.true)
  gr.tp <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.matrix) == 1)))
  gr.fp <- length(intersect(which(as.vector(cap.true.matrix) == 0), which(as.vector(cap.matrix) == 1)))
  gr.fn <- length(intersect(which(as.vector(cap.true.matrix) == 1), which(as.vector(cap.matrix) == 0)))
  gr.tn <- n*(n-1)/2-gr.tp-gr.fp-gr.fn
  gr.tpr <- gr.tp/(gr.tp+gr.fn)
  gr.fpr <- gr.fp/(gr.fp+gr.tn)
  gr.ri <- (gr.tp+gr.tn)/(n*(n-1)/2)
  gr.df <- data.frame('gr.tpr' = gr.tpr, 'gr.fpr' = gr.fpr, 'gr.ri' = gr.ri)
  return(gr.df)
}

spaperFun <- function(gamma, gamma.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate variable selection performance (RI, TPR, FPR)
  ## -----------------------------------------------------------------------------------------------------------------
  spa.tp <- length(intersect(which(as.vector(gamma.true) != 0), which(as.vector(gamma) != 0)))
  spa.tn <- length(intersect(which(as.vector(gamma.true) == 0), which(as.vector(gamma) == 0)))
  spa.fp <- length(intersect(which(as.vector(gamma.true) == 0), which(as.vector(gamma) != 0)))
  spa.tpr <- spa.tp/length(which(as.vector(gamma.true) != 0))
  spa.fpr <- spa.fp/length(which(as.vector(gamma.true) == 0))
  spa.ri <- (spa.tp+spa.tn)/length(gamma.true)
  spa.df <- data.frame('spa.tpr' = spa.tpr, 'spa.fpr' = spa.fpr, 'spa.ri' = spa.ri)
  return(spa.df)
}

mseFun <- function(beta, beta.true, gamma, gamma.true){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate mse
  ## -----------------------------------------------------------------------------------------------------------------
  mse.beta <- sqrt(sum((beta.true-beta)^2)/length(beta.true))
  mse.gamma <- sqrt(sum((gamma.true-gamma)^2)/length(gamma.true))
  mse.df <- data.frame('mse.beta' = mse.beta, 'mse.gamma' = mse.gamma)
  return(mse.df)
}

perFun <- function(result, wholedata){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Evaluate performance 
  ## -----------------------------------------------------------------------------------------------------------------
  group.true <- wholedata$group.true
  beta.true <- wholedata$beta.true
  gamma.true <- wholedata$gamma.true
  
  mse.beta <- mseFun(result$beta.gf, beta.true, result$gamma.gf, gamma.true)$mse.beta
  mse.gamma <- mseFun(result$beta.gf, beta.true, result$gamma.gf, gamma.true)$mse.gamma
  
  spa.per <- spaperFun(result$gamma.gf, gamma.true)
  spa.tpr <- spa.per$spa.tpr
  spa.fpr <- spa.per$spa.fpr
  spa.ri <- spa.per$spa.ri
  
  gr.num <- result$gr.num
  gr.gf <- result$gr.gf
  gr.per <- grperFun(gr.gf, group.true)
  gr.tpr <- gr.per$gr.tpr
  gr.fpr <- gr.per$gr.fpr
  gr.ri <- gr.per$gr.ri
  
  BIC.var <- result$BIC.var
  
  per.df <- data.frame('BIC.var' = BIC.var, 'gr.num' = gr.num, 'gr.ri' = gr.ri, 'spa.tpr' = spa.tpr, 'spa.fpr' = spa.fpr, 'spa.ri' = spa.ri, 
                       'gr.tpr' = gr.tpr, 'gr.fpr' = gr.fpr, 'mse.beta' = mse.beta, 'mse.gamma' = mse.gamma)
  
  return(list(per.df = per.df, gr.gf = gr.gf))
}

tuninglambda <- function(wholedata, xi.init, gamma.init, para, 
                         iter.max = 50, epsi = 0.01, mcp.para = 3, penal.para = 1, multi = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Tuning lambda parameters
  ## -----------------------------------------------------------------------------------------------------------------
  lambda1 <- para$lambda1; L1 <- length(lambda1)
  lambda2 <- para$lambda2; L2 <- length(lambda2)
  admmres <- vector(mode = "list", length = L1*L2)
  
  per.df <- data.frame()
  for (i in 1:L2) {
    for (j in 1:L1) {
      cat('-----------', i, '-th lambda2', lambda2[i] , 'and', j, '-th lambda1', lambda1[j] ,'--------------\n')
      admmres[[(i-1)*L1+j]] <- ADMM0(wholedata = wholedata, xi.init = xi.init, gamma.init = gamma.init, 
                                     lambda1 = c(0,lambda1[j]), lambda2 = c(0,lambda2[i]), tau = 0,
                                     xi.prior = xi.init, gamma.prior = gamma.init, admmstep = 2, 
                                     iter.max = iter.max, epsi = epsi, mcp.para = mcp.para, penal.para = penal.para, multi = multi)
      peri <- data.frame('tau' = 0, 'lambda1' = lambda1[j], 'lambda2' = lambda2[i], (perFun(admmres[[(i-1)*L1+j]], wholedata))$per.df)
      per.df <- rbind(per.df, peri)
    }
  }
  
  res.tune <- admmres[[min(which(per.df$BIC.var == min(per.df$BIC.var)))]]
  per.df.tune <- per.df[min(which(per.df$BIC.var == min(per.df$BIC.var))),]
  
  return(list(admmres = admmres, per.df = per.df, res.tune = res.tune, per.df.tune = per.df.tune))
}

tuningtau <- function(wholedata, xi.init, gamma.init, xi.prior, gamma.prior, sp, tau, lambda1.opt, lambda2.opt,
                      iter.max = 50, epsi = 0.01, mcp.para = 3, penal.para = 1, multi = F){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Tuning tau parameter with fixed lambda
  ## -----------------------------------------------------------------------------------------------------------------
  L0 <- length(tau)
  admmres <- vector(mode = "list", length = L0)
  
  per.df <- data.frame()
  for (i in 1:L0) {
    cat('-----------', i, '-th tau', tau[i] ,'--------------\n')
    admmres[[i]] <- ADMM0(wholedata = wholedata, xi.init = xi.init, gamma.init = gamma.init,
                          lambda1 = c(lambda1.opt,lambda1.opt), lambda2 = c(lambda2.opt,lambda2.opt), sp = sp, tau = tau[i],
                          xi.prior = xi.prior, gamma.prior = gamma.prior, admmstep = 2,
                          iter.max = iter.max, epsi = epsi, mcp.para = mcp.para, penal.para = penal.para, multi = multi)
    
    peri <- data.frame('tau' = tau[i], 'lambda1' = lambda1.opt, 'lambda2' = lambda2.opt, (perFun(admmres[[i]], wholedata))$per.df)
    per.df <- rbind(per.df, peri)
  }
  tau.opt <- tau[min(which(per.df$BIC.var == min(per.df$BIC.var)))]
  
  res.tune <- admmres[[min(which(per.df$BIC.var == min(per.df$BIC.var)))]]
  per.df.tune <- per.df[min(which(per.df$BIC.var == min(per.df$BIC.var))),]
  
  return(list(admmres = admmres, per.df = per.df, res.tune = res.tune, per.df.tune = per.df.tune))
}


###########################################################################################################################
## Part III ###############################################################################################################
##################################### Functions for calculating the initial values ########################################
initFun <- function(wholedata, lambda.sparse, div.num = 4, lambda1 = 1, lambda2 = 0.5){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Generate the initial value for proposed method
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  y.mat <- wholedata$data.y; n <- nrow(y.mat)
  x.mat <- wholedata$data.x; q <- ncol(x.mat)
  z.mat <- wholedata$data.z; p <- ncol(z.mat)
  x.row.vec <- vector(mode = 'list', length = n)
  z.row.vec <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    x.row.vec[[i]] <- matrix(x.mat[i,], nrow = 1, ncol = q)
  }
  for (i in 1:n) {
    z.row.vec[[i]] <- matrix(z.mat[i,], nrow = 1, ncol = p)
  }
  group.rough <- wholedata$group.rough
  l.matrix <- LMatrixFun(q, group.rough)
  a.matrix <- AMatrixFun(p, group.rough)
  x.mat.diag <- bdiag(x.row.vec)
  x.mat.diag.l <- x.mat.diag%*%l.matrix
  z.mat.diag <- bdiag(z.row.vec)
  
  # ridge
  q.matrix <- diag(n)-x.mat.diag.l%*%solve(t(x.mat.diag.l)%*%x.mat.diag.l)%*%t(x.mat.diag.l)
  gamma.ridge <- solve(t(z.mat.diag)%*%q.matrix%*%z.mat.diag+lambda1*diag(n*p)+lambda2*(t(a.matrix)%*%a.matrix))%*%t(z.mat.diag)%*%q.matrix%*%y.mat
  xi.ridge <- as.matrix(solve(t(x.mat.diag.l)%*%x.mat.diag.l)%*%t(x.mat.diag.l)%*%(y.mat-z.mat.diag%*%gamma.ridge))
  
  # compute xi.init
  xi.init <- as.matrix(xi.ridge)
  xi.init.mat <- matrix(xi.init, nrow = q)
  
  # compute refined groups
  gamma.ridge.mat <- matrix(gamma.ridge, nrow = p)
  gamma.ridge.mean <- apply(gamma.ridge.mat[1:7,], 2, mean)
  group.fine <- NULL
  for (k in 1:length(group.rough)) {
    kgroup <- group.rough[[k]]
    kn <- length(kgroup)
    kgamma.ridge.mean <- gamma.ridge.mean[kgroup]
    kgamma.order <- order(kgamma.ridge.mean)
    size <- kn/div.num
    for (i in 1:div.num) {
      ind <- kgamma.order[c(((i-1)*size+1):(i*size))]
      group.fine <- c(group.fine,list(sort(kgroup[ind])))
    }
  }
  
  # compute gamma.init
  y.mat <- as.matrix(y.mat-x.mat.diag%*%LMatrixFun(q, group.rough)%*%xi.init)
  mBIC <- NULL
  gamma <- NULL
  for (i in 1:length(lambda.sparse)) {
    alpha.list <- NULL
    for (j in 1:length(group.fine)) {
      jgroup <- group.fine[[j]]
      jy.mat <- as.matrix(y.mat[jgroup,], ncol = 1)
      jz.mat <- as.matrix(z.mat[jgroup,], ncol = p)
      lm.lasso <- glmnet(x = jz.mat, y = jy.mat, family = "gaussian", alpha = 1, penalty.factor = rep(1,p))
      alpha.list <- c(alpha.list, list(unname(coef(lm.lasso, s = lambda.sparse[i])[-1])))
    }
    alpha <- unlist(alpha.list)
    gammai <- as.matrix(LMatrixFun(p, group.fine)%*%alpha)
    gamma <- cbind(gamma, gammai)
    mBIC <- c(mBIC, log(sum((y.mat-as.matrix(z.mat.diag)%*%gammai)^2)/n)+log(n*p)*log(n)*length(which(alpha != 0))/n)
  }
  ind <- min(which(mBIC == min(mBIC)))
  lambda.sparse.select <- lambda.sparse[ind]
  gamma.init <- gamma[,ind]
  gamma.init.mat <- matrix(gamma.init, nrow = p, ncol = n)
  
  return(list(xi.init = xi.init, gamma.init = gamma.init, 
              xi.init.mat = xi.init.mat, gamma.init.mat = gamma.init.mat,
              group.fine = group.fine, lambda.sparse = lambda.sparse.select))
}


