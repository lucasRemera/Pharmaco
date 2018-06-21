###############
# medicaments #
###############

######################
# expected frequency #
######################
library(reshape2)
expectedMatrix=function(m){ #where m is the contingence matrix
  N=sum(m)
  meltedcontingence=expand.grid(atc=apply(m, 1, sum),dysp=apply(m, 2, sum))
  #M0=matrix(meltedcontingence[,2],nrow=nrow(contingence))
  M0=matrix(apply(meltedcontingence, 1, prod)/(N),nrow=nrow(m))
  return(M0)
}
####################
## Gibbs sampling ##
####################

O=10 #observed cases
E=1 #expected cases

modelM<-"
model{
A1~dgamma(alpha1,beta1)
A2~dgamma(alpha2,beta2)
lambda<-P*A1+(1-P)*A2
tau<-lambda*E
O~dpois(tau)

alpha1~dexp(.1)
beta1~dexp(0.1)
alpha2~dexp(.1)
beta2~dexp(.1)
P~dbeta(5,95)
}
"
library(rjags)
require(coda)
modelM.spec<-textConnection(modelM)
jags <- jags.model(modelM.spec,
                   data = list('O' = O,
                               'E' = E),
                   inits = list(alpha1=.2,beta1=.1,alpha2=2,beta2=4,P=0.5),
                   n.chains=4, 
                   n.adapt=100)
update(jags,1000)
samps.coda <- coda.samples(jags,
                           c('lambda','P'),
                           n.iter=100000
)
plot((samps.coda[[1]][,c("P")]))
quantile((samps.coda[[1]][,c("lambda")]),probs = c(.025,.5,.975))

###################################
## exact posterior distribution  ##
###################################
dq=function(x,alpha=1,beta=1,E=1){
  #print((beta))
  #return((1+(beta/E))**(-x)*(1+E/beta)**(-alpha)*gamma(alpha+x)/(gamma(alpha)*factorial(x))  )
  A1=(1+(beta/E))**(-x)*(1+E/beta)**(-alpha)
  A2=lgamma(alpha+x)-lgamma(alpha)-lfactorial(x)
  return(A1*exp(A2))
} #the density of probabilty of f=Pr(N=n), dumouchel 1999

pq=function(q,alpha=1,beta=1,E=1){
  if(length(q)==1) return(sum(sapply(0:q, function(i) dq(i,alpha,beta,E))))
  else return(sapply(q,function(qq) pq(qq,alpha,beta,E)))
} #the repartition function of f=Pr(N=n), dumouchel 1999

rq=function(n,alpha=1,beta=1,E=1){
  if(n==1){
    t=runif(1,0,1)
    nn=0
    q=pq(nn,alpha,beta,E)
    while(q<t){
      ##nn2=nn+1
      ##q2=pq(nn2,alpha,beta,E)
      #if(q2<t) nn=nn2
      #else if((q+q2)/2<t) nn=nn2
      ##q=q2
      
      nn=nn+1
      q=pq(nn,alpha,beta,E)
    }
    return(nn)
  }
  else return(replicate(n,{rq(1,alpha,beta,E)}))
} #simulate variables following the 'dq' law


rpi=function(n,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  P*rgamma(n,alpha1,beta1)+(1-P)*rgamma(n,alpha2,beta2)
} #simulate variables following the 'pi' law, dumouchel

rqn=function(n,E,O,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  #A1=0
  #A2=0
  #while(A1==0&A2==0){
  A1=rq(n,alpha1,beta1,E)
  A2=rq(n,alpha2,beta2,E)
  #print(paste0(A1,"  ",A2))
  qq = P*A1/(P*A1+(1-P)*A2)
  #}
  idx=which(A1==0&A2==0)
  while(length(idx)>0){
    AA1=rq(length(idx),alpha1,beta1,E)
    AA2=rq(length(idx),alpha2,beta2,E)
    A1[idx]=AA1
    A2[idx]=AA2
    idx=which(A1==0&A2==0)
  }
  qq = P*A1/(P*A1+(1-P)*A2)
  return(qq)
} #simulate variables following the 'Qn' law, dumouchel

a1=.2
b1=.1
a2=2
b2=4
pp=1/3 #the prior distribution from dumouchel 1999
E=10 #just an exemple of 'O'bserved and 'E'xpected data
O=100
rlambda=function(n,E,O,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  vqn=rqn(n,E,O,alpha1,beta1,alpha2,beta2,P)
  #print(vqn)
  return(rpi(n,alpha1+O,beta1+E,alpha2+O,beta2+E,vqn))
} #simulate lambda, the risk ratio

vlambda=rlambda(1000,E,O,a1,b1,a2,b2,pp)