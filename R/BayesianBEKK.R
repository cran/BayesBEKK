BayesianBEKK<-function(X,sd=10,iter=100,burnIn=40)
 {
    pop<-X
  k=dim(pop)[2]
  nT=dim(pop)[1]
  RTN<-pop[,1:k]
  COV1=cov(RTN)
  m1=BEKK11(pop)
  par<-m1$estimates
  AO=matrix(c(par[2],par[3],0,par[4]),2,2)
    AOAOt=AO%*%t(AO)
  A1=matrix(c(par[6:9]),2,2)
    B1=matrix(c(par[10:13]),2,2)
    resi=cbind(RTN[,1]-par[1],RTN[,2]-par[2])

  Sig=matrix(c(COV1),1,4)

  res=resi%*%t(A1)

  A1at<-function(x)
  {
    x=as.matrix(x)
    Prod1=NULL
    for(i in 1:nrow(x))
    {
      Prod1=rbind(Prod1,kronecker(x[i,],x[i,]))
    }
    Prod1
  }
  ArchP=A1at(res)

  llike=0
  for(t in 2:nT)
  {
    Sigt=AOAOt+matrix(ArchP[t-1,],k,k)+B1%*%
      matrix(Sig[t-1,],k,k)%*%t(B1)
    Sigt=(Sigt+t(Sigt))/2
    Sig=rbind(Sig,c(Sigt))
    dl=dmvnorm(resi[t,],mean=rep(0,k),Sigt,log=TRUE)
    llike=llike-dl
  }
  llike

  prior <- function(param)
  {
    mu1 = param[1]
    mu2 = param[2]
    a011 = param[3]
    a021 = param[4]
    a022 = param[5]
    a11 = param[6]
    a21 = param[7]
    a12 = param[8]
    a22 = param[9]
    b11 = param[10]
    b21 = param[11]
    b12 = param[12]
    b22 = param[13]

    mu1prior = dnorm(mu1, sd, log = T)
    mu2prior = dnorm(mu2, sd, log = T)
    a011prior = dnorm(a011, sd, log = T)
    a021prior = dnorm(a021, sd, log = T)
    a022prior = dnorm(a022, sd, log = T)
    a11prior = dnorm(a11, sd, log = T)
    a21prior = dnorm(a21, sd, log = T)
    a12prior = dnorm(a12, sd, log = T)
    a22prior = dnorm(a22, sd, log = T)
    b11prior = dnorm(b11, sd, log = T)
    b21prior = dnorm(b21, sd, log = T)
    b12prior = dnorm(b12, sd, log = T)
    b22prior = dnorm(b22, sd, log = T)


    return(mu1prior + mu2prior + a011prior+a021prior + a022prior + a11prior +a21prior +
             a12prior + a22prior + b11prior + b21prior + b12prior + b22prior)
  }

  posterior <- function(param)
  {
    return (llike + prior(param))
  }

  proposalfunction <- function(param){
    return(rnorm(13,mean = param,sd))
  }

  run_metropolis_MCMC <- function(startvalue, iterations){
    chain = array(dim = c(iterations+1,13))
    chain[1,] = startvalue
    for (i in 1:iterations){
      proposal = proposalfunction(chain[i,])

      probab = exp(posterior(proposal) - posterior(chain[i,]))
      if (runif(1) < probab){
        chain[i+1,] = proposal
      }else{
        chain[i+1,] = chain[i,]
      }
    }
    return(mcmc(chain))
  }
  chain= run_metropolis_MCMC(par, iter)
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
  return(summary(chain))
}


