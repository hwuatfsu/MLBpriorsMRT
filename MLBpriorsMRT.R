library(MfUSampler)
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
  
  beta.gibbs = matrix(0,P,B)
  lambda.gibbs = matrix(1,P,B)
  tau.gibbs = matrix(1,1,B)
  tauB.gibbs = matrix(1,1,B)
  alpha.gibbs = matrix(1,1,B)
  kappa.gibbs = matrix(2,1,B)
  alpha.beta.gibbs = matrix(1,1,B)
  kappa.beta.gibbs = matrix(2,1,B)
  
  for (b in 2:B){
    
    
    #update beta
    alpha.update = rbind(alpha.gibbs[b-1]*matrix(1,n1,1), y2+0.5 ,(r+0.5)*matrix(1,n3,1), alpha.beta.gibbs[b-1]*matrix(1,P,1))
    kappa.update = rbind(kappa.gibbs[b-1]*matrix(1,n1,1), (N+1)*matrix(1,n2,1), (y3+r+1), kappa.beta.gibbs[b-1]*matrix(1,P,1))
    H.update = rbind(-tau.gibbs[b-1]*X1star,X2star,X3star,tauB.gibbs[b-1]*diag(lambda.gibbs[,b-1]))
    mu.B.update = rbind(-tau.gibbs[b-1]*y1,matrix(0,(n2+n3+P),1))
    betasim = rbeta((n+P),alpha.update,kappa.update-alpha.update)
    w = log(betasim/(1-betasim))
    beta.gibbs[,b] = solve(t(H.update)%*% H.update,tol = 1e-22)%*%t(H.update)%*% mu.B.update + solve(t(H.update)%*% H.update, tol = 1e-22)%*%t(H.update)%*%w
    
    #update tau
    u_tau <- (tau.gibbs[b-1]^(a.tau-1+n1))*exp(-rexp(1,1))
    v_tau <- ((1+tau.gibbs[b-1])^(-a.tau-b.tau))*exp(-rexp(1,1))
    
    lower_tau <- u_tau^(1/(a.tau-1+n1))
    upper_tau <- v_tau^(-1/(a.tau+b.tau))-1
    
    
    if (tau.gibbs[b-1] >lower_tau ){
      ini_tau <- tau.gibbs[b-1]
    } else {
      ini_tau <- lower_tau
    }
    log.tau <- function(tau,beta,alpha,kappa){
      return(alpha*sum(tau*(y1-X1star%*%beta))-
               kappa*sum(log(matrix(1,n1,1)+exp(tau*(y1-X1star%*%beta)))))
      
    }
    tau.gibbs[b] <- MfU.Sample(x = ini_tau, f = log.tau, uni.sampler = "slice", alpha=alpha.gibbs[b-1],
                               kappa = kappa.gibbs[b-1], beta = beta.gibbs[,b],
                               control = MfU.Control(n=1,slice.lower = lower_tau,slice.upper = upper_tau))
    
    #update alpha.gibbs and kappa.gibbs
    alpha.gibbs[b] <- alpha.gibbs[b-1]
    
    tryCatch({
      f<-function(x){
        return(-n1*(lgamma(x))-n1*(lgamma(kappa.gibbs[b-1]-x))+x*sum(tau.gibbs[b]*(y1-X1star%*%beta.gibbs[,b])))
      }
      fprima<-function(x){
        return(-n1*digamma(x)+n1*digamma(kappa.gibbs[b-1]-x)+sum(tau.gibbs[b]*(y1-X1star%*%beta.gibbs[,b])))
      }
      alphastar<-seq(0.001,kappa.gibbs[b-1]-0.001,by = (kappa.gibbs[b-1])/(1000-1))
      likelihoods<- f(alphastar)
      ind<-which(likelihoods==max(likelihoods))
      simalpha<-ars2(1,f,fprima,x=c(0.99*alphastar[ind],alphastar[ind],kappa.gibbs[b-1]-0.001),lb=TRUE,xlb=1,ub=TRUE,xub =kappa.gibbs[b-1])
    }, error=function(e){})
    
    if (simalpha>0){
      alpha.gibbs[b] <- simalpha
    }
    kappa.gibbs[b] <- alpha.gibbs[b] *2
    tryCatch({
      f<-function(x){
        return(n1*lgamma(x)-n1*lgamma(x-alpha.gibbs[b])-x*sum(log(1+exp(tau.gibbs[b]*(y1-X1star%*%beta.gibbs[,b])))))
      }
      fprima<-function(x){
        return(n1*digamma(x)-n1*digamma(x-alpha.gibbs[b])-sum(log(1+exp(tau.gibbs[b]*(y1-X1star%*%beta.gibbs[,b])))))}
      kappastar=seq(alpha.gibbs[b]+1,4000,by = (4000-alpha.gibbs[b]-1)/(1000-1))
      likelihoods = f(kappastar)
      ind=which(likelihoods==max(likelihoods))
      kappa.gibbs[b]<-ars2(1,f,fprima,x=c(alpha.gibbs[b]+1,kappastar[ind],kappastar[ind]+100),lb=TRUE,xlb=alpha.gibbs[b]+1)
    }, error=function(e){})
    
    
    log.lam <- function(lambda,beta,alpha.beta,kappa.beta,tauB){
      
      return(alpha.beta*sum(tauB*lambda*beta)-
               kappa.beta*sum(log(1+exp(tauB*lambda*beta))))
      
    }
    #update lambda
    for (j in 1:P){
      u_lam <- (lambda.gibbs[j,b-1]^a.lam)*exp(-rexp(1,1))
      v_lam <- ((1+lambda.gibbs[j,b-1])^(-a.lam-b.lam))*exp(-rexp(1,1))
      
      lower <- min(u_lam^(1/a.lam),1)
      upper <- min(v_lam^(-1/(a.lam+b.lam))-1,10)
      
      if (lambda.gibbs[j,b-1] >lower ){
        ini_lam <-lambda.gibbs[j,b-1]
      } else {
        ini_lam <- lower
      }
      #library(crch)
      lambda.gibbs[j,b]<- MfU.Sample(x = ini_lam, f = log.lam, uni.sampler = "slice",
                                     beta = beta.gibbs[j,b],tauB = tauB.gibbs[,b-1],
                                     alpha.beta = alpha.beta.gibbs[b-1],kappa.beta =  kappa.beta.gibbs[b-1],
                                     control = MfU.Control(n=1,slice.lower = lower,slice.upper = upper))
    }
    
    
    
    #update tauB
    u_tauB <- (tauB.gibbs[b-1]^(a.taub+P-1))*exp(-rexp(1,1))
    v_tauB <- ((1+tauB.gibbs[b-1])^(-a.taub-b.taub))*exp(-rexp(1,1))
    
    lower_tauB <- u_tauB^(1/(a.taub+P-1))
    upper_tauB <- v_tauB^(-1/(a.taub+b.taub))-1
    
    if (tauB.gibbs[b-1] >lower_tauB ){
      ini_tauB <- tauB.gibbs[b-1]
    } else {
      ini_tauB <- lower_tauB
    }
    log.tauB <- function(tauB,beta,alpha.beta, kappa.beta,lambda){
      return(alpha.beta*sum(tauB*diag(lambda)%*%beta)-
               kappa.beta*sum(log(matrix(1,P,1)+exp(tauB*diag(lambda)%*%beta))))
      
    }
    tauB.gibbs[b] <- MfU.Sample(x = ini_tauB, f = log.tauB, uni.sampler = "slice",
                                beta = beta.gibbs[,b], alpha.beta=alpha.beta.gibbs[b-1],
                                kappa.beta=  kappa.beta.gibbs[b-1],lambda=lambda.gibbs[,b],
                                control = MfU.Control(n=1,slice.lower = lower_tauB,slice.upper = upper_tauB))
    
    
    #update alpha.beta and kappa.beta
    alpha.beta.gibbs[b]=alpha.beta.gibbs[b-1]
    
    tryCatch({
      f<-function(x){
        return(-P*lgamma(x)-P*lgamma(kappa.beta.gibbs[b-1]-x)+x*sum(tauB.gibbs[b]*diag(lambda.gibbs[,b])*beta.gibbs[,b]))}
      
      fprima<-function(x){
        return(-P*digamma(x)+P*digamma(kappa.beta.gibbs[b-1]-x)+sum(tauB.gibbs[b]*diag(lambda.gibbs[,b])*beta.gibbs[,b]))}
      alphabetastar=seq(0.001,kappa.beta.gibbs[b-1]-1,by = (kappa.beta.gibbs[b-1]-1)/(1000-1))
      likelihoods = f(alphabetastar)
      ind=which(likelihoods==max(likelihoods))
      simalphb<-ars2(1,f,fprima,x=c(0.99*alphabetastar[ind],alphabetastar[ind],kappa.beta.gibbs[b-1]-1),lb=TRUE,xlb=1,ub=TRUE,xub =kappa.beta.gibbs[b-1])
    }, error=function(e){})
    
    if (simalphb>0){
      alpha.beta.gibbs[b]= simalphb
    }
    kappa.beta.gibbs[b]= alpha.beta.gibbs[b]*2
    tryCatch({
      f<-function(x){
        return(P*lgamma(x)-P*lgamma(x-alpha.beta.gibbs[b])-x*sum(log(1+exp(tauB.gibbs[b]*diag(lambda.gibbs[,b])*beta.gibbs[,b]))))}
      fprima<-function(x){
        return(P*digamma(x)-P*digamma(x-alpha.beta.gibbs[b])-sum(log(1+exp(tauB.gibbs[b]*diag(lambda.gibbs[,b])*beta.gibbs[,b]))))}
      kappabetastar=seq(alpha.beta.gibbs[b]+1,4000,by = (4000-alpha.beta.gibbs[b]-1)/(1000-1))
      likelihoods = f(kappabetastar)
      ind=which(likelihoods==max(likelihoods))
      simkappb<-ars2(1,f,fprima,x=c(alpha.beta.gibbs[b]+1,kappabetastar[ind],kappabetastar[ind]+100),lb=TRUE,xlb=alpha.beta.gibbs[b]+1)
    }, error=function(e){})
    boundary.check=0
    if (simkappb>alpha.beta.gibbs[b]){
      boundary.check=1
    }
    
    if (boundary.check==1){
      kappa.beta.gibbs[b] <- simkappb
    }
    
    
    
    if (b%%1000 == 0)
    {
      print(b)
    }
    
  }
  
  betahat <-rowMeans(beta.gibbs[,burnin:B])
  #RMSE
  rmse_mlb_all <- sqrt(t((Xstar%*%betastar) - Xstar%*%betahat)%*%(Xstar%*%betastar - Xstar%*%betahat)/n)

  shared_th_signal3 <- MLBselect(beta.gibbs,B=5000, burnin = 2000,th=0.32, method = "intervals", alpha = 0.05)

  fpr_mlb <- fpr(shared_th_signal3,betastar)
  fnr_mlb <- fnr(shared_th_signal3,betastar)
  

  
  list <- list(rmse_mlb_all, fpr_mlb, fnr_mlb)
  colnames(list) <- c("rmse_mlb_all",  "fpr_hgt", "fnr_hgt")
  
  return(list)
  
}

