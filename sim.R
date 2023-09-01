# This code provides the simulation study from Wu and Bradley (2023).
# The model assumes that data y are of multiple response-type, where y1 is continuous-valued, y2 is integer-valued, and y3 is count-valued.
# X1star, X2star, and X3star are design matrices, such that y has a has a linear mixed effects representation on the appropriate link scale denoted with Xstar*beta.
# B The number of iterations of the block Gibbs sampler with burnin.
# N The sample size associated with the binomial observations.
# r The sample size associated with the negative binomial observations.
# a.lam and b.lam The values for the shape parameter associated with lambda.
# a.tau and b.tau The values for the shape parameter associated with tau.
# a.taub and b.taub The values for the shape parameter associated with tau_beta.
# The function returns root mean squared error (RMSE), average false negative rate (FNR), and average false positive rate (FPR).

library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)

comp_sim <- function(seed, p){
  
  n <- 100
  n1 <- 35
  n2 <- 35
  n3 <- 30
  index1 <- 1:n1
  index2 <- (n1+1):(n1+n2)
  index3 <- (n1+n2+1):(n)
  P <- 3*p
  
  set.seed(100)
  Xt <- matrix(runif(n*p),n,p)
  X <- Xt
  X[,1] <- sin(pi*Xt[,1]*Xt[,2])
  X[,2] <- (Xt[,3] - 0.5)^2
  X <- scale(X, center = TRUE, scale = TRUE)
  X1 <- X[index1,]
  X2 <- X[index2,]
  X3 <- X[index3,]
  
  X.add <-  as.matrix(t(apply(X[,c(4,5,6,7)], 1, combn, 2, prod)))
  R <- ncol(X.add)
  X.add1 <- X.add[index1,]
  X.add2 <- X.add[index2,]
  X.add3 <- X.add[index3,]
  
  X.add1.12 <- rbind(X.add1[18:35,], X.add2[1:17,])
  X.add1.13 <- rbind(X.add1[1:17,], X.add3[13:30,])
  
  X.add2.12 <- rbind(X.add1[18:35,], X.add2[1:17,])
  X.add2.23 <- rbind(X.add2[18:35,], X.add3[1:17,])
  
  X.add3.13 <- rbind(X.add1[1:15,], X.add3[16:30,])
  X.add3.23 <- rbind(X.add2[21:35,], X.add3[1:15,])
  
  X1star <- cbind(X1, matrix(0, n1, p), matrix(0, n1, p), X.add1.12, X.add1.13, matrix(0, n1, R))
  X2star <- cbind(matrix(0, n2, p), X2, matrix(0, n2, p), X.add2.12, matrix(0, n2, R), X.add2.23)
  X3star <- cbind(matrix(0, n3, p),matrix(0, n3, p), X3,matrix(0, n3, R),X.add3.13, X.add3.23)
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
  
  set.seed(seed)
  sigma <- 1
  N <- 40
  r <- 20
  y1 <- as.matrix(rnorm(n1, w1star, sigma))
  y2 <- as.matrix(rbinom(n2, size=N, prob=exp(w2star)/(1+exp(w2star))))
  y3 <- as.matrix(rnbinom(n3, size=r, prob=exp(w3star)/(1+exp(w3star))))
  y <- as.matrix(c(y1,y2,y3))
  
  
  B = 5000
  burnin = 2000
  a.lam =0.5
  b.lam= 0.5
  
  a.tau=10
  b.tau =10
  
  a.taub =10
  b.taub = 10
  

  out <- MLBpriorsMRT(y1,y2,y3,X1star,X2star,X3star,B,burnin,a.lam,b.lam,a.tau,b.tau,a.taub,b.taub,N,r)
  
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



registerDoParallel(8) 

re10 <- foreach (seed=100:119,.combine=rbind) %dopar% {
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


sim.p <- as.factor(c(rep(10,b),rep(20,b), rep(30,b), rep(40,b)))
sim.rmse <- as.numeric(c(re10[,1], re20[,1], re30[,1], re40[,1]))
sim.fpr <- as.numeric(c(re10[,2], re20[,2], re30[,2], re40[,2]))
sim.fnr <- as.numeric(c(re10[,3], re20[,3], re30[,3], re40[,3]))
sim.dat <- data.frame(sim.p = sim.p, sim.rmse=sim.rmse, sim.fpr=sim.fpr, sim.fnr=sim.fnr)


boxrmse <- ggplot(sim.dat, aes(x=sim.p, y=sim.rmse)) + 
  geom_boxplot(outlier.shape = NA)+
  coord_cartesian()+
  theme(text = element_text(size=20),
        axis.text.x = element_text(angle=90, hjust=1))+
  xlab("p") + ylab("RMSE")

boxfnr<-ggplot(sim.dat, aes(x=sim.p, y=sim.fnr)) + 
  geom_boxplot(outlier.shape = NA)+
  theme(text = element_text(size=25),
        axis.text.x = element_text(angle=90, hjust=1))+
  xlab("p") + ylab("FNR")

ggarrange(boxrmse,boxfnr, 
          ncol = 2, nrow = 1)
