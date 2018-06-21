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

dqn=function(x,E,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  #A1=0
  #A2=0
  #while(A1==0&A2==0){
  A1=dq(x,alpha1,beta1,E)
  A2=dq(x,alpha2,beta2,E)
  #print(paste0(A1,"  ",A2))
  qq = P*A1/(P*A1+(1-P)*A2)
  return(qq)
} #the 'Qn' law (density), dumouchel

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

# a1=.2
# b1=.1
# a2=2
# b2=4
# pp=1/3 #the prior distribution from dumouchel 1999
# E=10 #just an exemple of 'O'bserved and 'E'xpected data
# O=100
rlambda=function(n,E,O,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  vqn=rqn(n,E,O,alpha1,beta1,alpha2,beta2,P)
  #print(vqn)
  return(rpi(n,alpha1+O,beta1+E,alpha2+O,beta2+E,vqn))
} #simulate lambda, the risk ratio

# vlambda=rlambda(1000,E,O,a1,b1,a2,b2,pp)
# quantile(vlambda,probs = c(0.025,0.5,0.975)) #median + 95% credibility interval

################
# EB statistic #
################

EB=function(O,E,alpha1,beta1,alpha2,beta2,P){
  if(length(O)==1){
    qq=dqn(O,E,alpha1,beta1,alpha2,beta2,P)
    eb=qq*(digamma(alpha1+O)-log(beta1+E))+(1-qq)*(digamma(alpha2+O)-log(beta2+E))
    return(eb)
  }else{
    return(sapply(1:length(O), function(i) EB(O[i],E[i],alpha1,beta1,alpha2,beta2,P )))
  }
  
}


LLtheta=function(X,o=O,e=E,minus=-1){
  a1=X[1]
  b1=X[2]
  a2=X[3]
  b2=X[4]
  pp=X[5]
  nn=length(o)
  if(pp<0 | pp>1 | a1<=0 | a2<=0 | b1<=0 | b2<=0) ll=Inf
  else ll=sum(log(pp)+log(1-pp)+log(dq(o,a1,b1,e))+log(dq(o,a2,b2,e)))
  return(minus*ll)
}

mO=contingence #tableau de contingence des associations observees
mE=expectedMatrix(mO)
mmO=melt(mO)
mmE=melt(mE)

X0=c(.2,.05,.5,.5,.0004)
opt=optim(X0,LLtheta,o=mmO$value,e=mmE$value) #estimation by likelihood maximum
X1=opt$par


## par expectation maximization
LLtheta2=function(X,o=O,e=E,pp=.5,minus=-1,vpi){
  a1=X[1]
  b1=X[2]
  a2=X[3]
  b2=X[4]
  nn=length(o)
  if(pp<0 | pp>1 | a1<=0 | a2<=0 | b1<=0 | b2<=0) ll=Inf
  else ll=sum(vpi*log(dq(o,a1,b1,e)*pp)+(1-vpi)*log(dq(o,a2,b2,e)*(1-pp)))
  return(minus*ll)
}
PP=.1
XX=X0[1:4]
nsim=10
for(i in 1:nsim){
  print(i)
  pp1=sapply(1:nrow(mmO), function(ii) log(PP)+log(dq(mmO$value[ii], XX[1],XX[2],mmE$value[ii] )) )
  pp2=sapply(1:nrow(mmO), function(ii) log(1-PP)+log(dq(mmO$value[ii], XX[3],XX[4],mmE$value[ii] )) )
  ppi=exp(pp1)/(exp(pp1)+exp(pp2))
  #ppi=ifelse(pp1>pp2,1,0)
  #ppi=exp(pp1)
  PP=mean(ppi)
  XX=optim(XX,LLtheta2,o=mmO$value,e=mmE$value,pp=PP,vpi=ppi)$par
  
}

mm[which(ppi>.5),]
X1=c(XX,PP)

#E=1 #just an exemple of 'O'bserved and 'E'xpected data
#O=10
#vlambda=rlambda(1000,E,O,X1[1],X1[2],X1[3],X1[4],X1[5])
#hist(vlambda)
ebpost=EB(mmO$value,mmE$value,X1[1],X1[2],X1[3],X1[4],X1[5])
mm=cbind(mmO,mmE$value,ebpost)
mm[order(mm$ebpost,decreasing = T),] #ranking des associations

vlambda=rlambda(1000,mm[24662,4],mm[24662,3],X1[1],X1[2],X1[3],X1[4],X1[5])
quantile(vlambda,probs=c(.025,.5,.975)) #test if the association is significant (Monte-Carlo)
