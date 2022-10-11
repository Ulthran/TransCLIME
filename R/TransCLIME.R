#library(mvtnorm)
library(fastclime)
#library(Matrix)
library(glasso) # For joint graph estimator (Guo et al 2011)
library(lavaSearch2) # Only to symmetrize matrices
library(caret)
#library(BDcocolasso) # admm for spd projection

#' Algorithm based on fastclime package
#' min \|Theta\|_1
#' subject to \|cov(X)%*%Theta-Bmat\|_max <=lambda (if X is the raw sample)
#' subject to \|X%*%Theta-Bmat\|_max <=lambda (if X is the sample covariance matrix)
#'
#' @param X is SOMETHING
#' @param Bmat is SOMETHING
#' @param lambda is SOMETHING (default: 0.1)
#' @param scale is SOMETHING (default: T)
#' @param n is SOMETHING
#'
#' @return is a named list containing \itemize{
#' \item \code{Theta.hat} SOMETHING
#' \item \code{conv} SOMETHING
#' }
#'
#' @importFrom stats cov
#' @importFrom stats cov2cor
#' @importFrom utils capture.output
Myfastclime.s<-function(X,Bmat,lambda=0.1, scale=T, n){
  p<-ncol(X)
  obj=rep(-1,2*p)
  obj_bar=rep(0,2*p)
  rhs_bar<-rep(1, 2*p)
  if(isSymmetric(X, tol=10^(-4))){
    Sig.hat<-X
  }else{
    Sig.hat<-cov(X)
  }
  Sig.hat0<-Sig.hat
  feasible=T
  Theta.hat<-NULL
  Sig.diag<-diag(Sig.hat)
  Sig.hat<-cov2cor(Sig.hat)
  mat=rbind(cbind(Sig.hat,-Sig.hat),
              cbind(-Sig.hat,Sig.hat))
  for(j in 1:p){
    rhs <- c(Bmat[,j],-Bmat[,j])
    out.txt<-capture.output(  fastlp.re<-fastclime::fastlp(obj=obj, mat=mat, rhs=rhs+rhs_bar*lambda))
    if(!grepl("optimal", out.txt) ){
        feasible=F
        break
    }
    Theta.hat<-cbind(Theta.hat,(fastlp.re[1:p]-fastlp.re[-(1:p)])/sqrt(Sig.diag[j])/sqrt(Sig.diag))
    if(scale & sum(Theta.hat[,j]==0)==p){
      feasible=F
      break
    }else if(scale){
      Theta.hat[,j]<-  as.numeric(Bmat[j,j]/ (Sig.hat0[j,]%*%Theta.hat[,j]))*Theta.hat[,j]
     # Theta.hat[,j]<-as.numeric(Theta.hat[j,j]/ (t(Theta.hat[,j])%*%Sig.hat0%*%Theta.hat[,j]))*Theta.hat[,j]

    }
  }
  if(!feasible){
    cat('Theta.hat not found','\n')
    Theta.hat<-solve(cov(X)+diag(lambda,p))%*%Bmat
  }
  list(Theta.hat=Theta.hat, conv=feasible)
}

#' Compute the semi-positive definite projection of SigA.hat
#' with smallest eigenvalue lower bounded by eps
#'
#' @param SigA.hat is the input matrix
#' @param eps is the epsilon value
#'
#' @return is a named list containing \itemize{
#' \item \code{mat} SOMETHING
#' \item \code{conv} SOMETHING
#' }
Spd.proj<-function(SigA.hat, eps=NULL){
  p=ncol(SigA.hat)
  if(is.null(eps)){
    eps<-5/p
  }
  feasible=1
  SigA.t<-SigA.hat
   if(min(eigen(SigA.t)$values) <=eps ){
     feasible=2
     SigA.t<-ADMM_proj(SigA.hat, epsilon=eps)$mat
   }
  SigA.t<-lavaSearch2:::symmetrize(SigA.t, update.upper = TRUE)
  list(mat=SigA.t, conv=feasible)
}

#' Sourced from: https://rdrr.io/github/celiaescribe/BDcocolasso/src/R/ADMM_proj.R
#' @param v is SOMETHING
#' @param b is SOMETHING
#'
#' @return is SOMETHING
l1proj<-function(v, b){

  # Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
  # Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML

  stopifnot(b>0)

  u <- sort(abs(v),decreasing=TRUE)
  sv <- cumsum(u)
  rho <- max(which(u>(sv-b)/1:length(u)))
  theta <- max(0, (sv[rho]-b)/rho)
  w <-sign(v) * pmax(abs(v)-theta,0)

  return(w)
}

#' Sourced from: https://rdrr.io/github/celiaescribe/BDcocolasso/src/R/ADMM_proj.R
#' ADMM algorithm
#'
#' Finds the nearest positive semi-definite matrix with respect to the max norm
#'
#' @param mat Matrix to be projected
#' @param epsilon Approximation of the space of positive semi-definite matrix at epsilon
#' @param mu Penalty parameter of ADMM algorithm
#' @param it.max Number maximum of iterations
#' @param etol Tolerance parameter for the convergence of primal and dual residual
#' @param etol_distance Tolerance parameter for the convergence of the distance
#'
#' @return list containing \itemize{
#' \item \code{mat} Projected matrix
#' \item \code{df_ADMM} dataframe containing parameters of the ADMM algorithm for each iteration of the algorithm
#' }
#'
#' @seealso \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
ADMM_proj<-function(mat,
                    epsilon=1e-4,
                    mu=10,
                    it.max=1e3,
                    etol=1e-4,
                    etol_distance = 1e-4){



  p<-nrow(mat)

  # Initialization
  R<-diag(mat)
  S<-matrix(0,p,p)
  L<-matrix(0,p,p)

  itr<-0
  iteration <- eps_R <- eps_S <- eps_primal <- time <- distance <- NULL
  while (itr<it.max) {
    #print(itr)
    Rp<-R
    Sp<-S
    start <- Sys.time()
    # Subproblem I: R step
    W<-mat+S+mu*L
    W.eigdec<-eigen(W, symmetric=TRUE)
    W.V<-W.eigdec$vectors
    W.D<-W.eigdec$values
    R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)

    # Subproblem II: S step
    M<-R-mat-mu*L
    S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)
    for (i in 2:p){
      for (j in 1:(i-1)){
        S[j,i]<-S[i,j]
      }
    }

    # L step: update the Lagrange parameter
    L<-L-(R-S-mat)/mu
    end <- Sys.time()
    #Stocking the values of different parameters with the number of iterations
    iteration <- c(iteration, itr)
    eps_R <- c(eps_R,max(abs(R-Rp)))
    eps_S <- c(eps_S,max(abs(S-Sp)))
    eps_primal <- c(eps_primal, max(abs(R-S-mat)))
    time <- c(time, end - start)
    distance <- c(distance,max(abs(R-mat)))

    # Stopping Rule
    #cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if (((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)) || (abs(max(abs(Rp-mat)) - max(abs(R-mat)))<etol_distance)){
      itr<-it.max
    } else {
      itr<-itr+1
    }

    if (itr%%20==0) {
      mu<-mu/2
    }
  }
  df_ADMM <- data.frame(iteration = iteration, eps_R = eps_R, eps_S=eps_S, eps_primal=eps_primal, time=time, distance=distance)
  return(list(mat=R,df_ADMM=df_ADMM))

}

#' The main TransCLIME algorithm
#'
#' @param X is the primary data set
#' @param X.A is the set of X^{(k)} for k in A
#' @param const is SOMETHING
#' @param agg is a boolean to perform LS aggregation or not (default: T)
#' @param X.til is the sample set for aggregation (default: NULL)
#' @param Theta.cl is the CLIME estimator
#'
#' @return is SOMETHING
#'
#' @export
#'
#' @importFrom stats cov
#' @importFrom stats sd
Trans.CLIME<-function(X,X.A, const, agg=T, X.til=NULL,Theta.cl){
  if(agg &is.null(X.til)){
    cat('no aggregation samples provided.','\n')
  }
  n0=nrow(X)
  nA<-nrow(X.A)
  p<-ncol(X)
  sigA.hat<-mean(apply(X.A, 2, sd))
  sig0.hat<-mean(apply(X, 2, sd))

  lam.delta<-2*sig0.hat*sqrt(log(p)/n0)
  omega.l1<-mean(apply(Theta.cl,2, function(x) sum(abs(x))))
  Delta.re <- Myfastclime.s(X=diag(1,p), Bmat=diag(1,p)-t(Theta.cl)%*%cov(X.A), lambda=omega.l1*sqrt(log(p)/n0) , scale=F)
  if(Delta.re$conv){Delta.hat<-Delta.re$Theta.hat}else{ Delta.hat<-diag(0,p) }
  Theta.re <- Myfastclime.s(X=cov(X.A), Bmat=diag(1,p)-t(Delta.hat),
                            lambda=2*const*sqrt(log(p)/nA))
  Theta.hat<-Theta.re$Theta.hat

  if(agg){
    Omega.hat<-Agg(Theta.init=cbind(Theta.cl, Theta.hat), X.til=X.til)
  }else{
    Omega.hat<-Theta.hat
  }
  Omega.hat
}

#' LS aggregation function
#'
#' @param Theta.init is a column-bound combination of Theta.cl and Theta.hat
#' @param X.til is a subset of primary samples
#'
#' @return is SOMETHING
#'
#' @importFrom stats cov
Agg<-function(Theta.init, X.til){
  p<-ncol(X.til)
  n.til<-nrow(X.til)
  v.mat<-sapply(1:p, function(j){
    W.j<-cov(X.til%*%cbind(Theta.init[,j], Theta.init[,p+j]))
    v0=rep(0,2)
    v0[which.min(c(W.j[1,1]-2*Theta.init[j,j], W.j[2,2]-2*Theta.init[j,p+j]))]<-1
    v0
  })
  Theta.hat<-sapply(1:p, function(j) cbind(Theta.init[,j], Theta.init[,p+j])%*% v.mat[,j])

  Theta.hat
}



#' Cross validation for selecting tuning parameters
#'
#' @param X is SOMETHING
#' @param nfold is the number of folds (default: 5)
#'
#' @return is SOMETHING
#'
#' @importFrom caret createFolds
cv.clime<-function(X, nfold=5){
  p<-ncol(X)
  folds<-caret::createFolds(1:nrow(X), k=nfold)
  te<-NULL
  lam.seq<-seq(0.3,1.2,length.out=10)*2*sqrt(log(p)/nrow(X)*nfold/(nfold-1))
  for(i in 1:nfold){
    te<-rbind(te,sapply(lam.seq, function(lam){
      cur.clime<-Myfastclime.s(X=X[-folds[[i]],], Bmat=diag(1,p), lambda=lam)$Theta.hat
      Dist(X.test=X[folds[[i]],], Theta.hat=cur.clime, diag(1,ncol(X)))$te})
    )
  }

  te.ave<-colMeans(te)
  te.min<-which.min(te.ave)
  cat(te.ave,'\n')
  lam<-seq(0.3,1.2,length.out=10)[te.min]

  lam
}

#' Compute the estimation errors based on the test samples
#'
#' @param X.test is the set of test samples
#' @param Theta.hat is SOMETHING
#' @param Theta0 is SOMETHING
#'
#' @return is a named list containing \itemize{
#' \item \code{Frob} is SOMETHING
#' \item \code{S} is SOMETHING
#' \item \code{te} is SOMETHING
#' }
#'
#' @importFrom stats cov
Dist<- function(X.test, Theta.hat,Theta0){ ###compute the ell_2 error
  p<-ncol(Theta.hat)
  Theta.hat<-lavaSearch2:::symmetrize(Theta.hat, update.upper = TRUE)
  Theta<-Spd.proj(SigA.hat=Theta.hat, eps=0.001)$mat
  eigens<-eigen(Theta)$values
  te=sum(diag(cov(X.test)%*%Theta))/2-sum(log(eigens[eigens>0]))/2
  list(Frob=sum((Theta.hat-Theta0)^2)/p, S=max(abs(svd(Theta.hat-Theta0)$d))^2, te=te)
}

#' DB.clime.FDR is SOMETHING
#'
#' @param Theta is SOMETHING
#' @param X is SOMETHING
#'
#' @return is SOMETHING
#'
#' @importFrom stats cov
#' @importFrom stats pnorm
DB.clime.FDR<- function(Theta, X){
  n<-nrow(X)
  p<-ncol(X)
  Sig.hat<-cov(X)
  pval.all<-NULL
  diag.est<-rep(0,p)
  for(i in 1:p){
    diag.est[i]<-as.numeric(2*Theta[i,i] - t(Theta[,i])%*%Sig.hat%*%Theta[,i])
  }
  for( i in 1: (p-1)){
    for(j in (i+1):p){
      db.est<-as.numeric(Theta[i,j] +Theta[j,i]- t(Theta[,i])%*%Sig.hat%*%Theta[,j])
      # db.sd<- sd((X%*%Theta[,i])*(X%*%Theta[,j]))/sqrt(n)
      #db.sd<-sqrt(diag.est[i]*diag.est[j]+db.est^2)/sqrt(n)
      db.sd<-sqrt(Theta[i,i]*Theta[j,j]+Theta[i,j]^2)/sqrt(n)
      pval.ij<-2*(1-pnorm(abs(db.est/db.sd)))
      pval.all<-rbind(pval.all, c(i,j,pval.ij, db.est/db.sd))
    }
  }
  pval.all
}

#' FDR control function with input
#' @param p.val0 is the absolute value of z-scores
#' @param alpha is the fdr level
#' @param plot is a boolean to plot the output or not (default: F)
#'
#' @return is an array of indices for which z.abs is gte t.hat
#'
#' @importFrom graphics abline
#' @importFrom stats pnorm
#' @importFrom stats qnorm
BH.func<-function(p.val0, alpha, plot=F){
  print("p.val0")
  M=length(p.val0)
  Z.w<-qnorm(1-p.val0/2)
  fdr.est<-NULL
  t.seq<-seq(0,sqrt(2*log(M)-2*log(log(M))),0.01)
  for(t in t.seq){
    fdr.est<-c(fdr.est,M*2*(1-pnorm(t))/max(sum(Z.w>= t),1))
  }
  t.hat<-NULL
  t.hat<- t.seq[which(fdr.est<=alpha)[1]]
  if(plot){
    plot(t.seq,fdr.est, type='l')
    abline(h=alpha, col='blue')
  }

  if(is.na(t.hat)){
    t.hat<-sqrt(2*log(M))
  }
  #cat(t.hat, which(Z.w > t.hat),'\n')
  which(Z.w >= t.hat)
}

#' Joint graph estimator function from Guo et al. (2011)
#'
#' @param X.all is SOMETHING
#' @param n.vec is SOMETHING
#' @param lam.const is SOMETHING (default: NULL)
#'
#' @return is a named list containing \itemize{
#' \item \code{Theta.hat} is SOMETHING
#' \item \code{lam.const} is SOMETHING
#' }
#'
#' @importFrom stats cov
#'
#' @seealso \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3412604/}
jgl.fun<-function(X.all,n.vec, lam.const=NULL){
  K=length(n.vec)
  p <- ncol(X.all) # ADDED: p didn't exist in this context, is this a good way to initialize it? No idea
  X.list<-list()
  X.list[[1]] <- X.all[1:n.vec[1],]
  for(k in 2:K){
    ind.k<-(sum(n.vec[1:(k-1)])+1):(sum(n.vec[1:k]))
    X.list[[k]]<-X.all[ind.k,]
  }

  #initialization
  Theta.init<-list()
  for(k in 1: K){
    nu=2*sqrt(log(p)/n.vec[k])
    Theta.init[[k]]<-solve(cov(X.list[[k]])+diag(nu,p))
  }
  add <- function(x) Reduce("+", x)
  Theta.abs<- lapply(Theta.init, abs)
  weight<- sqrt(add(Theta.abs))
  weight<-apply(weight,1:2, function(x) 1/max(x,10^(-10)))
  ###tuning parameter
  if(is.null(lam.const)){
    bic.re<-NULL
    for(const in seq(0.2,2,length.out=10)){
      bic.cur<-0
      Theta.hat<-list()
      for(k in 1: K){
        Theta.re<-glasso::glasso(cov(X.list[[k]]), rho=const*2*sqrt(log(p)/n.vec[k])*weight, wi.init=Theta.init[[k]], maxit=100)
        Theta.hat[[k]]<-Theta.re$wi
        bic.cur <- bic.cur+Dist(X.test=X.list[[k]], Theta.hat[[k]], diag(1,p))$te+log(n.vec[k])*sum(Theta.hat[[k]]!=0)/2/n.vec[k]
      }
      bic.re<-c(bic.re,bic.cur)
    }
    lam.const=seq(0.2,2,length.out=10)[which.min(bic.re)]
    cat(lam.const,'\n')
  }
  ####running joint minimization for 10 rounds
  Theta.hat<-list()
  for(tt in 1:10){
    for(k in 1: K){
      Theta.re<-glasso::glasso(cov(X.list[[k]]), rho=lam.const*2*sqrt(log(p)/n.vec[k])*weight, wi.init=Theta.init[[k]], maxit=100)
      Theta.hat[[k]]<-Theta.re$wi
    }
    cat('tt=',tt,max(abs(Theta.hat[[1]]-Theta.init[[1]])),'\n')
    if(max(abs(Theta.hat[[1]]-Theta.init[[1]]))<=0.01){ break
    }else{
      Theta.init<-Theta.hat
    }
    add <- function(x) Reduce("+", x)
    Theta.abs<- lapply(Theta.init, abs)
    weight<- sqrt(add(Theta.abs))
    weight<-apply(weight,1:2, function(x) 1/max(x,10^(-10)))
  }
  return(list(Theta.hat=Theta.hat[[1]], lam.const=lam.const))
}
