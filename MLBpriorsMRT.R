library(MfUSampler)

  
  MLBpriorsMRT<-function(B,burnin,a.lam,b.lam,a.tau,b.tau,a.taub,b.taub,n,P,r){
  
  
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

  betahat<-rowMeans(beta.gibbs[,burnin:B])
  output<-list(betahat=betahat,beta=beta.gibbs,lambda=lambda.gibbs,tau=tau.gibbs,tauB=tauB.gibbs,alpha=alpha.gibbs,kappa=kappa.gibbs,
               alphaB=alpha.beta.gibbs,kappaB=kappa.beta.gibbs) 
  return(output)
}


