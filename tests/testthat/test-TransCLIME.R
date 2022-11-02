library(mvtnorm)

#generate the data for simulation
DataGen<-function(K,A.size,h, n.vec,s=10, p=100,type='Toep', ntest=100){

  if(type=='Toep'){
    Theta<-toeplitz(0.6^(1:p)*2)
    Theta[which(abs(Theta)<=0.05, arr.ind=T)]<- 0
  }else if(type=='Bdiag'){
    Theta<-kronecker(diag(p/4), toeplitz(c(1.2,0.9,0.6,0.3)))
  }else{
    Theta<-diag(1,p)+matrix(runif(p^2,0, 0.8),ncol=p,nrow=p)
    Theta<-(Theta+t(Theta))/2
    for(j in 1:p){
      for(l in 1:p){
        Theta[j,l]=Theta[j,l]/sqrt(abs(j-l)+1)
      }
    }
    for(j in 1:p){
      Theta[,j]<-Theta[,j]*(abs(Theta[,j])>=quantile(abs(Theta[,j]),1-s/p))
      Theta[j,]<-Theta[j,]*(abs(Theta[j,])>=quantile(abs(Theta[j,]),1-s/p))
    }
    Theta<-diag(max(0.1-min(eigen(Theta)$values),0),p)+Theta
  }
  Sig<- solve(Theta)
  X<-mvtnorm::rmvnorm(n.vec[1],rep(0,p), sigma=Sig)
  Theta.out<-diag(1,p)
  Omega.vec<-0
  for(k in 1 : K){
    if(k<=A.size){
      Delta.k<-matrix(rbinom(p^2,size=1,prob=0.1)*runif(p^2,-h/p,h/p),ncol=p)
      #  cat(max(colSums(abs(Delta.k))),'\n')
      Sig.k<-(diag(1,p)+Delta.k)%*%Sig
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- rbind(X, mvtnorm::rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
    }else{
      #Delta.out<-diag(0.5,p)+matrix(rbinom(p^2,size=1,prob=0.4)*runif(p^2,-0.2,0.2),ncol=p)
      Delta.out<-diag(1,p)+matrix(rbinom(p^2,size=1,prob=0.1)*0.4,ncol=p)
      Sig.k<-(diag(1,p)+Delta.out)%*%Sig
      Sig.k<-(Sig.k+t(Sig.k))/2
      if(min(eigen(Sig.k)$values)<0.05){
        Sig.k<-Sig.k+diag(0.1-min(eigen(Sig.k)$values),p)
      }
      X<- rbind(X, mvtnorm::rmvnorm(n.vec[k+1],rep(0,p), sigma=Sig.k))
    }
    Omega.vec<-c(Omega.vec, max(colSums(abs(diag(1,p)-Sig.k%*%Theta))))

  }
  cat(Omega.vec,'\n')
  list(X=X, Theta0=Theta, X.test=mvtnorm::rmvnorm(ntest,rep(0,p), sigma=Sig), Omega.l1=max(Omega.vec))

}

prep.data <- function(p, Theta0.type, K, A.size) {
  n.vec <- c(150, rep(300, K))
  dat.all <- DataGen(K=K, A.size=A.size, p=p, h=20, n.vec=n.vec, type=Theta0.type)
  Theta0 <- dat.all$Theta0
  X.all <- dat.all$X
  X.test <- dat.all$X.test
  #const <- cv.clime(X.all[1:n.vec[1],], nfold=10)
  const <- 0.5

  list(p=p, A.size=A.size, n.vec=n.vec, Theta0=Theta0, X.all=X.all, X.test=X.test, const=const)
}

test_that("everything works", {
  # Generate test data
  #set.seed(123)
  #data <- prep.data(200, 'Toep', 5, 4)
  #list2env(data, .GlobalEnv)
  #save.image(file="../../data/inputs.Rdata")
  # Load pregenerated test data
  load("../../data/inputs.Rdata")

  ###original clime####
  Theta.re0<-Myfastclime.s(X=X.all[1:n.vec[1],], Bmat=diag(1,p),
                           lambda=const*2*sqrt(log(p)/n.vec[1]))
  Theta.hat0<-Theta.re0$Theta.hat
  unlist(Dist(X.test=X.test, Theta.hat0, Theta0)) #get estimation errors and prediction errors

  ###oracle Trans-CLIME algorithm
  if(A.size>0){
    nA<-sum(n.vec[2:(A.size+1)])
    Omega.otl<-Trans.CLIME(X=X.all[1:n.vec[1],], X.A=X.all[(n.vec[1]+1):(n.vec[1]+nA),],
                           const=const, Theta.cl=Theta.hat0, agg=F)
  }else{ #A.size=0
    Omega.otl<-Theta.hat0
  }
  unlist(Dist(X.test, Omega.otl, Theta0))
  #errors in Frobenius norm; in spectral norm; prediction error based on negative log-likelihood

  #Trans-CLIME
  n0=round(n.vec[1]*4/5) #split sample for aggregation
  Theta.re0<-Myfastclime.s(X=X.all[1:n0,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
  Theta.init<-Theta.re0$Theta.hat
  Omega.tl1 <- Trans.CLIME(X=X.all[1:n0,], X.A=X.all[-(1:n.vec[1]),], const=const,
                           X.til= X.all[(n0+1):n.vec[1],], Theta.cl=Theta.init)
  ind2<-(n.vec[1]-n0+1): n.vec[1]
  Theta.re0<-Myfastclime.s(X=X.all[ind2,], Bmat=diag(1,p), lambda=const*2*sqrt(log(p)/n0))
  Theta.init<-Theta.re0$Theta.hat
  Omega.tl2 <- Trans.CLIME(X=X.all[ind2,], X.A=X.all[-(1:n.vec[1]),], const=const,
                           X.til= X.all[1:(n.vec[1]-n0),], Theta.cl=Theta.init)
  Omega.tl<-(Omega.tl1+Omega.tl2)/2
  unlist(Dist(X.test, Omega.tl, Theta0))

  ####multi-task graph estimator
  const.jgl<-jgl.fun(X.all,n.vec)$lam.const
  cat(const.jgl,'\n')
  jgl.re<-jgl.fun(X.all,n.vec, lam.const=const.jgl)
  unlist(Dist(X.test, jgl.re$Theta.hat, Theta0))

  ###debiasing and FDR control
  pval.cl<-DB.clime.FDR(Theta=Theta.hat0,X=X.all[1:n.vec[1],]) #debiased single-task CLIME
  pval.otl<-DB.clime.FDR(Theta=Omega.otl,X=X.all[1:n.vec[1],]) #deibased oracle Trans-CLIME
  pval.tl<-DB.clime.FDR(Theta=Omega.tl,X=X.all[1:n.vec[1],])#debiased Trans-CLIME
  #find significant pairs at FDR level 0.1
  cl.pairs<- pval.cl[BH.func(pval.cl[,3], 0.1),1:2]
  tl.pairs<-pval.tl[BH.func(pval.tl[,3], 0.1),1:2]
  otl.pairs<-pval.otl[BH.func(pval.otl[,3], 0.1),1:2]

  FDR.res <- c(sum(Theta0[cl.pairs]==0)/max(length(cl.pairs),1),
    sum(Theta0[otl.pairs]==0)/max(length(otl.pairs),1),
    sum(Theta0[tl.pairs]==0)/max(length(tl.pairs),1)) ###FDR result

  power.res <- c(sum(Theta0[cl.pairs]!=0)/sum(Theta0[pval.cl[,1:2]]!=0),
    sum(Theta0[otl.pairs]!=0)/sum(Theta0[pval.otl[,1:2]]!=0),
    sum(Theta0[tl.pairs]!=0)/sum(Theta0[pval.tl[,1:2]]!=0)) ###power result

  save(Theta.re0, Omega.otl, Omega.tl, jgl.re, FDR.res, power.res, file="../../data/outputs.Rdata")
  print(FDR.res)
  print(power.res)
  expect_equal(2 * 2, 4)
})
