comp_sim <- function(p, seed){
  
  n <- 100
  n1 <- 35
  n2 <- 35
  n3 <- 30
  index1 <- 1:n1
  index2 <- (n1+1):(n1+n2)
  index3 <- (n1+n2+1):(n)
  P <- 3*p
  
  set.seed(seed)
  Xt <- matrix(runif(n*p),n,p)
  X <- Xt
  X[,1] <- sin(pi*Xt[,1]*Xt[,2])
  X[,2] <- (Xt[,3] - 0.5)^2
  X <- scale(X, center = TRUE, scale = TRUE)
  X1 <- X[index1,]
  X2 <- X[index2,]
  X3 <- X[index3,]
  
  X.add <-  as.matrix(t(apply(X[,c(4,5,6,7)], 1, combn, 2, prod)))
  r <- ncol(X.add)
  X.add1 <- X.add[index1,]
  X.add2 <- X.add[index2,]
  X.add3 <- X.add[index3,]
  X.add1.12 <- rbind(X.add1[18:35,], X.add2[1:17,])
  X.add1.13 <- rbind(X.add1[1:17,], X.add3[13:30,])
  
  X.add2.12 <- rbind(X.add1[18:35,], X.add2[1:17,])
  X.add2.23 <- rbind(X.add2[18:35,], X.add3[1:17,])
  
  X.add3.13 <- rbind(X.add1[1:15,], X.add3[16:30,])
  X.add3.23 <- rbind(X.add2[21:35,], X.add3[1:15,])
  
  X1star <- cbind(X1, matrix(0, n1, p), matrix(0, n1, p), X.add1.12, X.add1.13, matrix(0, n1, r))
  X2star <- cbind(matrix(0, n2, p), X2, matrix(0, n2, p), X.add2.12, matrix(0, n2, r), X.add2.23)
  X3star <- cbind(matrix(0, n3, p),matrix(0, n3, p), X3,matrix(0, n3, r),X.add3.13, X.add3.23)
  Xstar <- rbind(X1star, X2star, X3star)
  
  beta1 <- as.matrix(c(10,20,10,10,10,10,0*seq(1,(p-6))))
  beta2 <- as.matrix(c(0,1,1,0,0,1,0*seq(1,(p-6))))
  beta3 <- as.matrix(c(2,2,2,0,0,0,0*seq(1,(p-6))))
  
  beta1.add <- as.matrix(c(1,1,0,0,0,0))
  beta2.add <- as.matrix(c(0,1,1,0,0,0))
  beta3.add <- as.matrix(c(0,0,1,1,0,0))
  
  betastar <- rbind(beta1, beta2, beta3, beta1.add, beta2.add, beta3.add)
  
  w1star <- X1star%*%betastar
  w2star <- X2star%*%betastar
  w3star <- X3star%*%betastar
  wstar <- rbind(w1star, w2star, w3star)
  
  
  sigma <- 1
  N <- 40
  r <- 20
  y1 <- as.matrix(rnorm(n1, w1star, sigma))
  y2 <- as.matrix(rbinom(n2, size=N, prob=exp(w2star)/(1+exp(w2star))))
  y3 <- as.matrix(rnbinom(n3, size=r, prob=exp(w3star)/(1+exp(w3star))))
  y <- as.matrix(c(y1,y2,y3))
  
  
  B = 5000
  burnin = 2000
  a.lam =1
  b.lam= 2
  
  a.tau=1
  b.tau =2
  
  a.taub =1
  b.taub = 2
  
  n = length(y)
  P = ncol(Xstar)
  r = 40
  
  out <- MLBpriorsMRT(B,burnin,a.lam,b.lam,a.tau,b.tau,a.taub,b.taub,n,P,r)
  
  out$betahat
  betahat <-  out$betahat
  
  #RMSE
  rmse_mlb_all <- sqrt(t((Xstar%*%betastar) - Xstar%*%betahat)%*%(Xstar%*%betastar - Xstar%*%betahat)/n)
  
  #FPR FNR
  shared_th_signal <- MLBselect(out$beta,B=5000, burnin = 2000,th=0.32, method = "intervals", alpha = 0.05)
  fpr_mlb <- fpr(shared_th_signal,betastar)
  fnr_mlb <- fnr(shared_th_signal,betastar)
  
  
  list <- list(rmse_mlb_all, fpr_mlb, fnr_mlb)

  return(list)
  
}


library(foreach)
library(doParallel)
registerDoParallel(8) 

re10 <- foreach (seed=200:219,.combine=rbind) %dopar% {
  comp_sim(seed, p =10)
}

re20 <- foreach (seed=100:119,.combine=rbind) %dopar% {
  comp_sim(seed, p =20)
}

re30 <- foreach (seed=100:119,.combine=rbind) %dopar% {
  comp_sim(seed, p =30)
}
re40 <- foreach (seed=100:119,.combine=rbind) %dopar% {
  comp_sim(seed, p = 40)
}
