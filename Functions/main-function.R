#####################################################################################
## Fundamental supporting functions and ADMM function of proposed methods.
#####################################################################################

################################################################################################################################
## Part I ######################################################################################################################
############################################# Some fundamental supporting functions ############################################
grindFun <- function(i, group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Find group index of the i-th subject 
  ## -----------------------------------------------------------------------------------------------------------------
  for (k in 1:length(group)) {
    if(i %in% group[[k]]){
      grind <- k
    }
  }
  return(grind)
}

LMatrixFun <- function(q, group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Sparse matrix L in computation
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group)); K <- length(group)
  L <- Matrix(0, nrow = n, ncol = K, sparse = T)
  for (i in 1:n) {
    L[i,grindFun(i, group)] <- 1
  }
  L <- kronecker(L, bdiag(diag(q)))
  return(L)
}

AMatrixFun <- function(p, group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Sparse matrix A in computation
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group)); K <- length(group); D <- list()
  for (k in 1:K) {
    kgroup <- group[[k]]; m <- length(kgroup); D[[k]] <- list()
    if(m > 1){
      for (i in 1:(m-1)) {
        D[[k]][[i]] <- Matrix(0, nrow = m-i, ncol = n, sparse = T)
        D[[k]][[i]][,kgroup[i]] <- 1
        for (l in ((i+1):m)) {
          D[[k]][[i]][l-i,kgroup[l]] <- -1
        }
      }
    }
    D[[k]] <- do.call(rbind, D[[k]])
  }
  D <- do.call(rbind, D)
  A <- kronecker(D, bdiag(diag(p))) 
  return(A)
}

pairindFun <- function(group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate the pair indices of subjects by group structure
  ## -----------------------------------------------------------------------------------------------------------------
  P <- vector(mode = "list")
  if(length(group) == length(unlist(group))){
    P <- matrix(c(0,0), 2, 1)
  }
  else{
    for (k in 1:length(group)) {
      kgroup <- group[[k]]; m <- length(kgroup); P[[k]] <- vector(mode = "list")
      if(m > 1){
        for (i in 1:(m-1)) {
          P[[k]][[i]] <- rbind(rep(kgroup[i],(m-i)), kgroup[(i+1):m])
        }
      }
      P[[k]] <- do.call(cbind, P[[k]])
    }
    P <- do.call(cbind, P)
  }
  return(P)
}

pairmatFun <- function(group){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate index upper triangular matrix by group structure
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group)); P <- matrix(0, nrow = n, ncol = n)
  pairind <- pairindFun(group); sumdif <- ncol(pairind)
  if(length(group) == n){
    P <- P
  }
  else{
    for (j in 1:sumdif) {
      P[pairind[1,j],pairind[2,j]] <- 1
    }
  }
  return(P)
}

groupFun <- function(pairmat){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate a group structure (list) by the output of pairmatFun
  ## -----------------------------------------------------------------------------------------------------------------
  pairmat2 <- pairmat+t(pairmat); diag(pairmat2) <- 1
  group.unique.matrix <- unique(pairmat2); gr.num <- nrow(group.unique.matrix)
  group <- vector(mode = "list", length = gr.num)
  for (k in 1:gr.num) {
    for (i in 1:nrow(pairmat2)) {
      if(sum(abs(pairmat2[i,]-group.unique.matrix[k,])) == 0){
        group[[k]] <- c(group[[k]],i)
      }
    }
  }
  return(group)
} 

groupgfFun <- function(eta.gf, p, group.rough){
  ## -----------------------------------------------------------------------------------------------------------------
  ## Generate group structure by estimated coefficients
  ## -----------------------------------------------------------------------------------------------------------------
  n <- length(unlist(group.rough)); sumdif <- length(eta.gf)/p
  pairind <- pairindFun(group.rough)
  diff.gf <- apply(matrix(abs(eta.gf), nrow = p, ncol = sumdif), 2, sum)
  cap.gf.matrix <- matrix(0, nrow = n, ncol = n)
  if(length(which(diff.gf == 0)) == 0){
    gr.gf <- groupFun(cap.gf.matrix)
  }
  else if(length(which(diff.gf == 0)) == sumdif){
    gr.gf <- group.rough
  }
  else{
    pairind.gf <- as.matrix(pairind[,which(diff.gf == 0)])
    for(j in 1:ncol(pairind.gf)){
      cap.gf.matrix[pairind.gf[1,j],pairind.gf[2,j]] <- 1
    }
    gr.gf <- groupFun(cap.gf.matrix)
  }
  gr.num <- length(gr.gf)
  return(list(gr.gf = gr.gf, gr.num = gr.num))
}


###########################################################################################################################
## Part II ################################################################################################################
############################################ Functions for main algorithms ################################################
ADMM0 <- function(wholedata, xi.init, gamma.init, lambda1, lambda2, sp, 
                  tau, xi.prior, gamma.prior, admmstep,
                  iter.max, epsi, mcp.para, penal.para, multi = F){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The ADMM implementation of Step1 or Step2.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ wholedata: The input data analyzed.
  ## @ xi.init: The Initial values of regression coefficients of X.
  ## @ gamma.init: The Initial values of regression coefficients of Z.
  ## @ lambda1: A 2-vector, the tuning parameters controlling the number of sparse coefficients in Step1/2.
  ## @ lambda2: A 2-vector, the tuning parameter controlling the number of refined subgroup in Step1/2.
  ## @ sp: A vector of indices of important features via prior information.
  ## @ tau: A weight in Step2.
  ## @ xi.prior: Regression coefficients results by Step1.
  ## @ gamma.prior: Regression coefficients results by Step1.
  ## @ admmstep: int 1,2, means Step1/2 respectively.
  ## @ iter.max: int >= 2, maximum number of cycles of the ADMM algorithm, the default setting is 50.
  ## @ epsi: A float value, algorithm termination threshold.
  ## @ mcp.para: The regularization parameter in MCP, the default setting is 3.
  ## @ penal.para: The penalty parameter in ADMM algorithm, the default setting is 1.
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
  sumdif <- nrow(a.matrix)/p
  
  x.mat.diag <- bdiag(x.row.vec)
  x.mat.diag.l <- x.mat.diag%*%l.matrix
  z.mat.diag <- bdiag(z.row.vec)
  
  ## initialization
  delta.init <- gamma.init
  eta.init <- a.matrix%*%gamma.init
  v.init <- Matrix(0, nrow = p*n, ncol = 1, sparse = T) 
  u.init <- Matrix(0, nrow = p*sumdif, ncol = 1, sparse = T) 
  
  iter <- 1
  xi.est.list <- vector(mode = "list", length = iter.max); xi.est.list[[iter]] <- xi.init 
  gamma.est.list <- vector(mode = "list", length = iter.max); gamma.est.list[[iter]] <- gamma.init
  delta.est.list <- vector(mode = "list", length = iter.max); delta.est.list[[iter]] <- delta.init
  eta.est.list <- vector(mode = "list", length = iter.max); eta.est.list[[iter]] <- eta.init 
  v.est.list <- vector(mode = "list", length = iter.max); v.est.list[[iter]] <- v.init 
  u.est.list <- vector(mode = "list", length = iter.max); u.est.list[[iter]] <- u.init 
  rm(xi.init, gamma.init, delta.init, eta.init, v.init, u.init)
  gc()
  
  ## iteration
  while(iter <= iter.max){
    iter <- iter+1
    if(iter > iter.max){
      cat("The step iteration number exceed the maximum values!\n")
      cat("The eps is now", eps, "\n")
      break
    }
    
    # update gamma,xi
    if(admmstep == 2){
      y.mat <- (1-tau)*y.mat+tau*(x.mat.diag.l%*%xi.prior+z.mat.diag%*%gamma.prior)
    }
    if(!multi){
      q.matrix <- diag(n)-x.mat.diag.l%*%solve(t(x.mat.diag.l)%*%x.mat.diag.l)%*%t(x.mat.diag.l)
    }
    else{
      q.matrix <- diag(n)-x.mat.diag.l%*%solve(t(x.mat.diag.l)%*%x.mat.diag.l+0.001*diag(q*length(group.rough)))%*%t(x.mat.diag.l)
    }
    
    gamma.est.list[[iter]] <- solve(t(z.mat.diag)%*%q.matrix%*%z.mat.diag+penal.para*(diag(n*p)+t(a.matrix)%*%a.matrix))%*%
                                   (t(z.mat.diag)%*%q.matrix%*%y.mat+penal.para*(delta.est.list[[iter-1]]+t(a.matrix)%*%eta.est.list[[iter-1]])
                                    -v.est.list[[iter-1]]-t(a.matrix)%*%u.est.list[[iter-1]])
    
    if(!multi){
      xi.est.list[[iter]] <- solve(t(x.mat.diag.l)%*%x.mat.diag.l)%*%t(x.mat.diag.l)%*%(y.mat-z.mat.diag%*%gamma.est.list[[iter]])
    }
    else{
      xi.est.list[[iter]] <- solve(t(x.mat.diag.l)%*%x.mat.diag.l+0.001*diag(q*length(group.rough)))%*%t(x.mat.diag.l)%*%(y.mat-z.mat.diag%*%gamma.est.list[[iter]])
    }
    
    # mcp update delta
    delta.fix <- gamma.est.list[[iter]] + v.est.list[[iter-1]]/penal.para
    delta.fix <- matrix(delta.fix, nrow = p, ncol = n)
    delta.opt <- Matrix(0, nrow = p, ncol = n, sparse = T)
    if(admmstep == 1){
      for (i in 1:n) {
        for (j in 1:p) {
          if(j %in% setdiff(1:p, sp) & abs(delta.fix[j,i]) <= mcp.para*lambda1[1]){
            delta.opt[j,i] <- delta.fix[j,i]*max(0,1-(lambda1[1]/penal.para)/abs(delta.fix[j,i]))/(1-1/(mcp.para*penal.para))
          }
          else{delta.opt[j,i] <- delta.fix[j,i]}
        }
      }
    }
    if(admmstep == 2){
      for (i in 1:n) {
        for (j in 1:p) {
          if(abs(delta.fix[j,i]) <= mcp.para*lambda1[2]){
            delta.opt[j,i] <- delta.fix[j,i]*max(0,1-(lambda1[2]/penal.para)/abs(delta.fix[j,i]))/(1-1/(mcp.para*penal.para))
          }
          else{delta.opt[j,i] <- delta.fix[j,i]}
        }
      }
    }
    dim(delta.opt) <- c(p*n,1)
    delta.est.list[[iter]] <- delta.opt
    
    # mcp update eta
    eta.est.list[[iter]] <- Matrix(0, nrow = p*sumdif, ncol = 1, sparse = T)
    eta.fix1 <- (a.matrix%*%gamma.est.list[[iter]])+u.est.list[[iter-1]]/penal.para 
    one.matrix <- bdiag(rep(list(rep(1, p)), sumdif))
    eta.fix2 <- sqrt(t(one.matrix)%*%(eta.fix1^2)) 
    if(admmstep == 1){
      indleq <- which(eta.fix2 <= mcp.para*lambda2[1])
      indgeq <- which(eta.fix2 > mcp.para*lambda2[1])
      eta.fix3 <- as.matrix(apply(1-lambda2[1]/(penal.para*eta.fix2), 1, function(x) max(x,0))[indleq])
    }
    if(admmstep == 2){
      indleq <- which(eta.fix2 <= mcp.para*lambda2[2])
      indgeq <- which(eta.fix2 > mcp.para*lambda2[2])
      eta.fix3 <- as.matrix(apply(1-lambda2[2]/(penal.para*eta.fix2), 1, function(x) max(x,0))[indleq])
    }
    eta.fix4 <- as.vector(apply(eta.fix3, 1, function(x) rep(x,p)))
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
    
    # update v,u
    v.est.list[[iter]] <- v.est.list[[iter-1]]+penal.para*(gamma.est.list[[iter]]-delta.est.list[[iter]])
    u.est.list[[iter]] <- u.est.list[[iter-1]]+penal.para*(a.matrix%*%gamma.est.list[[iter]]-eta.est.list[[iter]])
    
    eps <- max(sqrt(sum((a.matrix%*%gamma.est.list[[iter]]-eta.est.list[[iter]])^2)/(sumdif*p)),
               sqrt(sum((gamma.est.list[[iter]]-delta.est.list[[iter]])^2))/(n*p))
    if(eps < epsi){
      cat("The iterations finish in", iter, "-th step.\n")
      break
    }
    
    xi.est.list[iter-1] <- list(NULL)
    gamma.est.list[iter-1] <- list(NULL)
    delta.est.list[iter-1] <- list(NULL)
    eta.est.list[iter-1] <- list(NULL)
    u.est.list[iter-1] <- list(NULL)
    v.est.list[iter-1] <- list(NULL)
  }
  
  ##################
  if(iter > iter.max){
    xi.gf <- as.matrix(xi.est.list[[iter-1]])
    beta.gf <- as.matrix(l.matrix%*%xi.gf)
    gamma.gf <- as.matrix(gamma.est.list[[iter-1]])
    delta.gf <- as.matrix(delta.est.list[[iter-1]])
    eta.gf <- as.matrix(eta.est.list[[iter-1]])
  }
  else{
    xi.gf <- as.matrix(xi.est.list[[iter]])
    beta.gf <- as.matrix(l.matrix%*%xi.gf)
    gamma.gf <- as.matrix(gamma.est.list[[iter]])
    delta.gf <- as.matrix(delta.est.list[[iter]])
    eta.gf <- as.matrix(eta.est.list[[iter]])
  }
  rm(xi.est.list, gamma.est.list, delta.est.list, eta.est.list, v.est.list, u.est.list)
  gc()
  
  gr.gf <- groupgfFun(eta.gf, p, group.rough)$gr.gf
  gr.num <- groupgfFun(eta.gf, p, group.rough)$gr.num
  
  delta.gf.mat <- matrix(delta.gf, nrow = p, ncol = n)
  alpha.gf.list <- vector(mode = 'list', length = gr.num)
  for (k in 1:gr.num) {
    alpha.gf.list[[k]] <- matrix(apply(matrix(delta.gf.mat[,gr.gf[[k]]], nrow = p), 1, mean), ncol = 1)
  }
  alpha.gf <- do.call(rbind,alpha.gf.list)
  alpha.gf[which(abs(alpha.gf) < 0.001)] <- 0
  gamma.gf <- as.matrix(LMatrixFun(p, gr.gf)%*%alpha.gf)
  
  s <- length(which(abs(alpha.gf) != 0))
  residual <- log(sum((y.mat-x.mat.diag%*%beta.gf-z.mat.diag%*%gamma.gf)^2)/n)
  BIC.var <- residual+log(n*(p+q))*log(n)*(length(group.rough)*q+s)/n
  
  return(list(xi.gf = xi.gf, beta.gf = beta.gf, alpha.gf = alpha.gf, gamma.gf = gamma.gf,
              eta.gf = eta.gf, BIC.var = BIC.var, gr.num = gr.num, gr.gf = gr.gf))
}

ADMM <- function(wholedata, xi.init, gamma.init, lambda1, lambda2, sp, tau, 
                 iter.max = c(50,50), epsi = c(0.01,0.01), mcp.para = 3, penal.para = 1, multi = F){
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The ADMM implementation of proposed method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## ADMM Step1
  res1 <- ADMM0(wholedata = wholedata, xi.init = xi.init, gamma.init = gamma.init, 
                lambda1 = lambda1, lambda2 = lambda2, sp = sp, admmstep = 1, 
                iter.max = iter.max[1], epsi = epsi[1], mcp.para = mcp.para, penal.para = penal.para, multi = multi)
  xi.prior <- res1$xi.gf
  gamma.prior <- res1$gamma.gf
  ## ADMM Step2
  res2 <- ADMM0(wholedata = wholedata, xi.init = xi.init, gamma.init = gamma.init, 
                lambda1 = lambda1, lambda2 = lambda2, sp = sp, tau = tau, 
                xi.prior = xi.prior, gamma.prior = gamma.prior, admmstep = 2, 
                iter.max = iter.max[2], epsi = epsi[2], mcp.para = mcp.para, penal.para = penal.para, multi = multi)
  return(res2)
}



