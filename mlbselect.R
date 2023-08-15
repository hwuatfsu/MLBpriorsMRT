# The function performs variable selection for the MLB priors.
MLBselect <- function(beta.gibbs,B, burnin, alpha, th, method){
  if (method=="intervals"){
    effsamp <- B - burnin
    left <- floor(alpha * effsamp/2)
    right <- ceiling((1 - alpha/2) * effsamp)
    BetaSort <- apply(beta.gibbs[,burnin:B], 1, sort, decreasing = F)
    left.points <- BetaSort[left, ]
    right.points <- BetaSort[right, ]
    signals <- as.numeric(1 - ((left.points <= 0) & (right.points >= 0)))
    return(signals)
  }
  
  if (method == "threshold"){
    betahat <-rowMeans(beta.gibbs[,burnin:B]) 
    signals <- as.numeric(betahat >= th)
    return(signals)
  }
}
#false positive rate
fpr <- function(signals,true){
  fpr <- sum(signals[true==0])/length(true[true==0])
  return(fpr)
}
#false negative rate
fnr <- function(signals, true){
  fnr <- sum(signals[true!=0]==0)/length(true[true!=0])
  return(fnr)
}
