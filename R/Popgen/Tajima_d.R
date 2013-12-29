mySFStest <- function(sfs1,divide,sfs2=NULL){
  ## INPUT: a vector of frequency spectrums, length = n-1
  ##        frequency class from 1 ~ n-1, counts
  ##        divide is a percentage (0,1) which groups the sites into two classes
  ## DO   : calculate various theta estimators and
  ##        do several tests based on these estimators, following ZK's work
  ## OUTPUT: to be determined

  ## calculate different estimators of theta
  n <- length(sfs1)+1
  n1 <- n-1
  S <- sum(sfs1)
  theta.w <- S / sum(1/1:n1)
  theta.pi <- sum(1:n1 * n1:1 * sfs1) / choose(n,2)
  theta1 <- sfs1[1]
  theta.h <- sum(1:n1 * 1:n1 * sfs1) / choose(n,2)

  ## Calculate Tajima's D
  d <- theta.pi-theta.w
  a1 = sum(1/1:n1); b1 = (n+1)/(3*(n-1)); c1=b1-1/a1; e1=c1/a1
  a2 = sum(1/(1:n1)^2); b2 = 2*(n^2+n+3)/(9*n*(n-1)); c2=b2-(n+2)/a1/n+a2/(a1^2); e2=c2/(a1^2+a2)
  vd.hat <- e1*S+e2*S*(S-1)
  D <- d/sqrt(vd.hat)
  
  ## 2. divide the counts into two categories according to "divide"
  j <- trunc(n*divide) 
  obs <- c(sum(sfs1[1:j]),sum(sfs1[-(1:j)])) # groups: 1:j and j:n1
  if(is.null(sfs2)){
    sfs0 <- 1/1:n1 ## sfs under neutral model
    exp.p <- c(sum(sfs0[1:j]),sum(sfs0[-(1:j)])) / sum(sfs0)
    if(sum(obs<5))
      chisq <- chisq.test(obs,p=exp.p,simulate=TRUE,B=10000)
    else
      chisq <- chisq.test(obs,p=exp.p,simulate=FALSE)
    ## 4. likelihood ratio
    L0=dbinom(obs[1],size=sum(obs),prob=exp.p[1])
    L1=dbinom(obs[1],size=sum(obs),prob=obs[1]/sum(obs))
    LR=L0/L1
  }
  else{
    obs2 <- c(sum(sfs2[1:j]),sum(sfs2[-(1:j)])) # groups: 1:j and j:n1
    if(sum(obs)<5 | sum(obs2)<5)
      chisq <- chisq.test(obs,obs2, simulate=TRUE, B=10000)
    else
      chisq <- chisq.test(obs,obs2,simulate=FALSE)
    ## 4. likelihood ratio
    L0=dbinom(obs[1],size=sum(obs),prob=obs2[1]/sum(obs2))
    L1=dbinom(obs[1],size=sum(obs),prob=obs[1]/sum(obs))
    LR=L0/L1
  }

  
  result <- list(table=
                 cbind(c(paste("<",divide,sep=""),
                         paste(">=",divide,sep="")),
                       paste(obs," (",round(obs/sum(obs),2),"%)",sep=""),
                       paste(chisq$exp," (",round(chisq$exp/sum(chisq$exp),2),"%)",sep="")),
                 chisq.stat=chisq$stat, chisq.pvalue=chisq$p.value,
                 Likelihood.rtio=LR,
                 Tajima.D=D)
  return(result)
}

mySFStest1 <- function(sfs1,neut=NULL){
  ## a simple categorical X^2 test
  ## 1. divide the polymorphism counts into three categories of roughtly equal size
  n1 <- length(sfs1)
  n <- n1+1
  if(is.null(neut)){
    neut <- 1/1:n1 ## neutral spectrum
  }
  cum.neut <- cumsum(neut)
  cum <- max(cum.neut)
  cut1 <- c(which.min(abs(cum.neut - cum/3)),
            which.min(abs(cum.neut - cum*2/3)))
  category <- as.factor(c(rep("low",cut1[1]),
                          rep("intermediate",cut1[2]-cut1[1]),
                          rep("high", n1-cut1[2])))
  obs <- by(sfs1,category,sum)
  exp <- by(neut,category,sum)
  if(sum(obs<5))
    simulate1=TRUE
  else
    simulate1=FALSE
  barplot(rbind(obs,exp/sum(exp)*sum(obs))[,3:1],beside=TRUE,legend.text=TRUE)
  return(chisq.test(obs,p=exp/sum(exp),simulate=simulate1))
}
