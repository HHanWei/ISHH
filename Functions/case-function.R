#####################################################################################
# Functions for real data analysis.
#####################################################################################
casepreFun <- function(x.mat, y.mat, lambda1, lambda2, iter.max = 50, epsi = 0.001, mcp.para = 3, penal.para = 1, seed = 1234){
  n <- nrow(y.mat); q <- ncol(x.mat)
  
  # compute beta init
  set.seed(seed)
  kmeansres <-  kmeans(x.mat, 3, nstart = 25)
  beta.init.mat <- matrix(0, nrow = q, ncol = n)
  for (k in 1:3) {
    ind <- which(unname(kmeansres$cluster) == k)
    kx.mat <- matrix(x.mat[ind,], ncol = q)
    ky.mat <- y.mat[ind,]
    kbeta <- unname(solve(t(kx.mat)%*%kx.mat+0.001*diag(q))%*%t(kx.mat)%*%ky.mat)
    beta.init.mat[,ind] <- matrix(rep(kbeta,length(ind)),nrow = q)
  }
  beta.init <- matrix(beta.init.mat, nrow = q*n, ncol = 1)
  
  # admm for beta
  a.matrix <- AMatrixFun(q,list(1:n))
  sumdif <- nrow(a.matrix)/q
  
  x.row.vec <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    x.row.vec[[i]] <- matrix(x.mat[i,], nrow = 1, ncol = q)
  }
  x.mat.diag <- bdiag(x.row.vec)
  
  eta.init <- a.matrix%*%beta.init
  u.init <- Matrix(0, nrow = q*sumdif, ncol = 1, sparse = T) 
  
  iter <- 1
  beta.est.list <- vector(mode = "list", length = iter.max); beta.est.list[[iter]] <- beta.init
  eta.est.list <- vector(mode = "list", length = iter.max); eta.est.list[[iter]] <- eta.init 
  u.est.list <- vector(mode = "list", length = iter.max); u.est.list[[iter]] <- u.init 
  rm(beta.init, eta.init, u.init)
  gc()
  
  while(iter <= iter.max){
    iter <- iter+1
    
    if(iter > iter.max){
      cat("The step iteration number exceed the maximum values!\n")
      cat("The eps is now", eps, "\n")
      break
    }
    
    beta.est.list[[iter]] <- solve(t(x.mat.diag)%*%x.mat.diag+penal.para*t(a.matrix)%*%a.matrix+lambda1*diag(n*q))%*%
      (t(x.mat.diag)%*%y.mat+penal.para*t(a.matrix)%*%eta.est.list[[iter-1]]-t(a.matrix)%*%u.est.list[[iter-1]])
    
    eta.est.list[[iter]] <- Matrix(0, nrow = q*sumdif, ncol = 1, sparse = T)
    eta.fix1 <- (a.matrix%*%beta.est.list[[iter]])+u.est.list[[iter-1]]/penal.para 
    one.matrix <- bdiag(rep(list(rep(1, q)), sumdif))
    eta.fix2 <- sqrt(t(one.matrix)%*%(eta.fix1^2)) 
    indleq <- which(eta.fix2 <= mcp.para*lambda2)
    indgeq <- which(eta.fix2 > mcp.para*lambda2)
    eta.fix3 <- as.matrix(apply(1-lambda2/(penal.para*eta.fix2), 1, function(x) max(x,0))[indleq])
    eta.fix4 <- as.vector(apply(eta.fix3, 1, function(x) rep(x,q)))
    if(length(indleq) < 2){
      indleq2 <- which(rowSums(as.matrix(one.matrix[,indleq])) != 0)
      indgeq2 <- which(rowSums(one.matrix[,indgeq]) != 0)
    }
    if(length(indgeq) < 2){
      indleq2 <- which(rowSums(one.matrix[,indleq]) != 0)
      indgeq2 <- which(rowSums(as.matrix(one.matrix[,indgeq])) != 0)
    }
    if(length(indleq) >= 2 & length(indgeq) >= 2){
      indleq2 <- which(rowSums(one.matrix[,indleq]) != 0)
      indgeq2 <- which(rowSums(one.matrix[,indgeq]) != 0)
    }
    eta.est.list[[iter]][indgeq2,] <- eta.fix1[indgeq2,]
    eta.est.list[[iter]][indleq2,] <- eta.fix4*eta.fix1[indleq2,]/(1-1/(mcp.para*penal.para))
    
    u.est.list[[iter]] <- u.est.list[[iter-1]]+penal.para*(a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])
    
    eps <- sqrt(sum((a.matrix%*%beta.est.list[[iter]]-eta.est.list[[iter]])^2)/(sumdif*q))
    if(eps < epsi){
      cat("The iterations finish in", iter, "-th step.\n")
      break
    }
    
    beta.est.list[iter-1] <- list(NULL)
    eta.est.list[iter-1] <- list(NULL)
    u.est.list[iter-1] <- list(NULL)
  }
  
  ##################
  if(iter > iter.max){
    beta.gf <- as.matrix(beta.est.list[[iter-1]])
    eta.gf <- as.matrix(eta.est.list[[iter-1]])
  }
  else{
    beta.gf <- as.matrix(beta.est.list[[iter]])
    eta.gf <- as.matrix(eta.est.list[[iter]])
  }
  rm(beta.est.list, eta.est.list, u.est.list)
  gc()
  
  gr.gf <- groupgfFun(eta.gf, q, list(1:n))$gr.gf
  gr.num <- groupgfFun(eta.gf, q, list(1:n))$gr.num
  
  residual <- log(sum((y.mat-x.mat.diag%*%beta.gf)^2)/n)
  BIC.var <- residual+log(n*q)*log(n)*gr.num*q/n
  
  beta.mat.gf <- matrix(beta.gf, nrow = q)
  xi.gf <- matrix(0, nrow = q, ncol = length(gr.gf))
  for (k in 1:gr.num) {
    xi.gf[,k] <- apply(matrix(beta.mat.gf[,gr.gf[[k]]],nrow = q), 1, mean)
  }
  
  return(list(group.rough = gr.gf, xi.gf = xi.gf, BIC.var = BIC.var))
}

caseinitFun <- function(wholedata, lambda.sparse, div.num = div.num, lambda1 = 1, lambda2 = 0.5){
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
  q.matrix <- diag(n)-x.mat.diag.l%*%solve(t(x.mat.diag.l)%*%x.mat.diag.l+lambda1*diag(length(group.rough)*q))%*%t(x.mat.diag.l)
  gamma.ridge <- solve(t(z.mat.diag)%*%q.matrix%*%z.mat.diag+lambda1*diag(n*p)+lambda2*(t(a.matrix)%*%a.matrix))%*%t(z.mat.diag)%*%q.matrix%*%y.mat
  xi.ridge <- as.matrix(solve(t(x.mat.diag.l)%*%x.mat.diag.l+lambda1*diag(length(group.rough)*q))%*%t(x.mat.diag.l)%*%(y.mat-z.mat.diag%*%gamma.ridge))
  
  # compute xi.init
  xi.init <- as.matrix(unname(xi.ridge))
  xi.init.mat <- matrix(xi.init, nrow = q)
  
  # compute refined groups
  gamma.ridge.mat <- matrix(gamma.ridge, nrow = p)
  gamma.ridge.med <- apply(gamma.ridge.mat, 2, median)
  group.fine <- NULL
  for (k in 1:length(group.rough)) {
    kgroup <- group.rough[[k]]
    kn <- length(kgroup)
    kgamma.ridge.med <- gamma.ridge.med[kgroup]
    kgamma.order <- order(kgamma.ridge.med)
    size <- floor(kn/div.num[k])
    for (i in 1:(div.num[k]-1)) {
      ind <- kgamma.order[c(((i-1)*size+1):(i*size))]
      group.fine <- c(group.fine,list(sort(kgroup[ind])))
    }
    group.fine <- c(group.fine,list(sort(kgroup[kgamma.order[c(((div.num[k]-1)*size+1):kn)]])))
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





