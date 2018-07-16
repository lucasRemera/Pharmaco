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

###########
# Poisson #
###########
poissonTest=function(o,e,rr=FALSE){
  if(length(o)==1){
    if(rr) return(c(o/e,1-ppois(o,e)))
    else return(1-ppois(o,e))
  } 
  else return(sapply(1:length(o), function(ii) poissonTest(o[ii],e[ii],rr)))
}


####################
## Gibbs sampling ##
####################

# O=10 #observed cases
# E=1 #expected cases
# 
# modelM<-"
# model{
# A1~dgamma(alpha1,beta1)
# A2~dgamma(alpha2,beta2)
# lambda<-P*A1+(1-P)*A2
# tau<-lambda*E
# O~dpois(tau)
# 
# alpha1~dexp(.1)
# beta1~dexp(0.1)
# alpha2~dexp(.1)
# beta2~dexp(.1)
# P~dbeta(5,95)
# }
# "
# library(rjags)
# require(coda)
# modelM.spec<-textConnection(modelM)
# jags <- jags.model(modelM.spec,
#                    data = list('O' = O,
#                                'E' = E),
#                    inits = list(alpha1=.2,beta1=.1,alpha2=2,beta2=4,P=0.5),
#                    n.chains=4, 
#                    n.adapt=100)
# update(jags,1000)
# samps.coda <- coda.samples(jags,
#                            c('lambda','P'),
#                            n.iter=100000
# )
# plot((samps.coda[[1]][,c("P")]))
# quantile((samps.coda[[1]][,c("lambda")]),probs = c(.025,.5,.975))

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

qn=function(x,E,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  #A1=0
  #A2=0
  #while(A1==0&A2==0){
  A1=dq(x,alpha1,beta1,E)
  A2=dq(x,alpha2,beta2,E)
  #print(paste0(A1,"  ",A2))
  qq = P*A1/(P*A1+(1-P)*A2)
  return(qq)
} #the 'Qn' law (density), dumouchel ##wrong !!!


dqn=function(x,E,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  if(length(x)==1){
    if(x==0) return(dq(0,alpha1,beta1,E))
    toIntegrate=function(t) dq(t,alpha1,beta1,E)*dq(P*t*(1/x-1)/(1-P))
    return(integrate(toIntegrate,0,Inf)[[1]])
  }
  else return(sapply(x, function(xx) dqn(xx,E,alpha1,beta1,alpha2,beta2,P )))
}


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

rlambda2=function(n,E,O,alpha1=1,beta1=1,alpha2=1,beta2=1,P=.5){
  qq=qn(O,E,alpha1,beta1,alpha2,beta2,P)
  #print(vqn)
  return(rpi(n,alpha1+O,beta1+E,alpha2+O,beta2+E,qq))
}
# vlambda=rlambda(1000,E,O,a1,b1,a2,b2,pp)
# quantile(vlambda,probs = c(0.025,0.5,0.975)) #median + 95% credibility interval

################
# EB statistic #
################

EB=function(O,E,alpha1,beta1,alpha2,beta2,P,quick.and.dirty=FALSE, toLog=!quick.and.dirty){
  if(length(O)==1){
    qq=dqn(O,E,alpha1,beta1,alpha2,beta2,P)
    eb=qq*(digamma(alpha1+O)-log(beta1+E))+(1-qq)*(digamma(alpha2+O)-log(beta2+E))/log(2)
    if(quick.and.dirty){
      if(toLog) return(log(2**eb*(exp(c(-2,0,2)/sqrt(1+O))))/log(2))
      else return(2**eb*(exp(c(-2,0,2)/sqrt(1+O))))
    } 
    else if(!toLog) return(2**eb)
    else return(eb)
  }else{
    return(sapply(1:length(O), function(i) EB(O[i],E[i],alpha1,beta1,alpha2,beta2,P,quick.and.dirty,toLog )) )
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

idx=which(mmO$value!=0)
# mmO=mmO[idx,]
# mmE=mmE[idx,]
#filtre sur les observes >0?
X0=c(.2,.05,.5,.5,.0004)
opt=optim(X0,LLtheta,o=mmO$value,e=mmE$value,control = list(maxit=1000)) #estimation by likelihood maximum
X1=opt$par

ebpost=EB(mmO$value,mmE$value,X1[1],X1[2],X1[3],X1[4],X1[5])
mm=cbind(mmO,mmE$value,ebpost)
mm$idx=1:nrow(mm)
mm[order(mm$ebpost,decreasing = T),] #ranking des associations

vlambda=rlambda2(10000,mm[145,4],mm[145,3],X1[1],X1[2],X1[3],X1[4],X1[5])
quantile(vlambda,probs=c(.025,.5,.975)) #test if the association is significant (Monte-Carlo)

###############################
##  expectation maximization ##
###############################
 
mm[which(ppi>.5),]
head(mm[order(ppi,decreasing = T),], sum(ppi>.5))
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


#################################################
## Parameters estimation with a gibbs sampling ##
#################################################
# here observed are name obs and expected expe
N=length(obs)
modelMAll<-"
model{
for(i in 1:N){
A1[i]~dgamma(alpha1,beta1)
A2[i]~dgamma(alpha2,beta2)
lambda[i]<-P*A1[i]+(1-P)*A2[i]
tau[i]<-lambda[i]*E[i]
O[i]~dpois(tau[i])
}
alpha1~dexp(.1)
beta1~dexp(0.1)
alpha2~dexp(.1)
beta2~dexp(.1)
P~dbeta(5,95)
}
"
library(rjags)
require(coda)
modelMall.spec<-textConnection(modelMAll)
jags <- jags.model(modelMall.spec,
                   data = list('O' = obs,
                               'E' = expe,
                               'N'=N),
                   inits = list(alpha1=.2,beta1=.1,alpha2=2,beta2=4,P=0.5),
                   n.chains=3, 
                   n.adapt=100)
update(jags,1000)
samps.coda <- coda.samples(jags,
                           c('lambda','P'),
                           n.iter=10000
)
plot((samps.coda[[1]][,c("lambda")]))
quantile((samps.coda[[1]][,c("lambda")]),probs = c(.025,.5,.975))

###############
# simulations #
###############
#k to modulate numbers of total observation (N'=k*N)
#n the number of simulation
simulatePoissonH0=function(expected,k=1,n=1){
  if(n==1) return(rpois(length(expected),k*expected))
  else replicate(n,{simulatePoissonH0(expected,k,1)})
}

simulatePoissonH1=function(expected,idx,rr=2,k=1,n=1){
  
  if(n==1){
    RR=rep(1,length(expected))
    RR[idx]=rr
    return(rpois(length(expected),k*expected*RR))
  } 
  else replicate(n,{simulatePoissonH1(expected,idx,rr,k,1)})
}

eH0=simulatePoissonH0(mmE$value,1,10)
ebH0=apply(eH0, 2, function(ee){
  X00=c(.2,.1,2,4,.1)
  opt0=optim(X00,LLtheta,o=ee,e=mmE$value) #estimation by likelihood maximum
  X10=opt0$par 
  
  return(EB(ee,mmE$value,X10[1],X10[2],X10[3],X10[4],X10[5]))
})


seqQ=seq(0,1,by=0.01)
# plot(quantile(ebpost,probs=seqQ),quantile(ebH0[,1],probs=seqQ),col=rainbow(length(seqQ)))
# abline(a = 0,b=1,col="blue")
ggplot()+geom_point(aes(x=quantile(ebpost,probs=seqQ),y=quantile(ebH0[,1],probs=seqQ),col=seqQ))+geom_abline(aes(slope=1,intercept=0))+coord_fixed()

mebH0=melt(ebH0)     
ggplot()+geom_histogram(data=mebH0,aes(x=value,group=Var2,fill=factor(Var2)),alpha=.5,bins=100)+geom_histogram(aes(x=ebpost),col="red",alpha=.5,bins=100)

txT=sapply(seq(-2,2,by=0.001), function(s) sum(ebpost<s) )/length(ebpost)
tx0=sapply(seq(-2,2,by=0.001), function(s) sum(ebH0[,1]<s) )/length(ebH0[,1])
ggplot()+geom_line(aes(x=txT,y=tx0),col=3+seq(-2,2,by=0.001))+geom_abline(aes(slope=1,intercept=0))


X100=apply(eH0, 2, function(ee){
  X00=c(.2,.1,2,4,.1)
  opt0=optim(X00,LLtheta,o=ee,e=mmE$value) #estimation by likelihood maximum
  X10=opt0$par
  return(X10)
})

##############################################
# lasso regression for confounding controls ##
##############################################
library(glmnet)
load("burtSIDPNG.RData",verbose = T)

burtATC=burtSIDPNG[,4:100]
burtDYSP=burtSIDPNG[,101:242]
burtDYSP=as.data.frame(apply(burtDYSP,2,as.numeric))
burtATC=(apply(burtATC,2,as.numeric))

lasso=glmnet(burtATC,as.factor(burtDYSP[,"Q60"]),family = "binomial",alpha=1,lambda = 10**seq(-10,6,by=.1))
plot(lasso,xvar="lambda")
lasso_cv=cv.glmnet(burtATC,burtDYSP[,"Q54"],alpha=1,lambda = 10**seq(-6,6,by=.1))
plot(lasso_cv)
lasso_cv$lambda[80]
lasso_cv$glmnet.fit$beta[,82]


##################################
## Bayesian logistic regression ##
##################################

load("burtSIDPNG.RData",verbose = T)
burtATC=burtSIDPNG[,4:100]
burtDYSP=burtSIDPNG[,101:242]
burtDYSP=apply(burtDYSP,2,as.numeric)
burtATC=apply(burtATC,2,as.numeric)


modellog<-"
model{
for(i in 1:Nobs){
dysp[i] ~ dbern(p[i])
p[i] <- 1 / (1 + exp(-z[i]))
z[i]<-sum(m[i,])
for(j in 1:Natc){
m[i,j]<-beta[j]*atc[i,j]

}
}
for(k in 1:Natc){
beta[k]~dnorm(0,tau[k])
tau[k]~dexp(gamma)
}
gamma~dexp(1)
}
"

Nobs=nrow(burtATC)
Natc=ncol(burtATC)
dysp=burtDYSP[,43]
modellog.spec<-textConnection(modellog)
jags <- jags.model(modellog.spec,
                   data = list('Nobs' = Nobs,
                               'Natc' = Natc,
                               'dysp'=dysp,
                               'atc'=burtATC),
                   inits = list(gamma=1),
                   n.chains=3, 
                   n.adapt=100)
update(jags,1000)
samps.codalog <- coda.samples(jags,
                              c('beta','gamma','tau'),
                              n.iter=10000
)

save(samps.codalog,file="GibbsLogit180712.RData")

#########################
# information component #
#########################

MutualInformation=function(o,e){
  N=sum(o)
  sum(log2(o/e)*o/N)
}

InformationComponent=function(o,e,ATC,MALF,atc,malf){
  i=which(ATC==atc&MALF==malf)
  log2(o[i]/e[i])
}

#MutualInformation(mmO$value,mmE$value)

#########
# BCPNN #
#########
#keeping the Bate,1998 notation
BCPNN=function(cxy,cx,cy,C,aalpha1=1,aalpha=2,bbeta1=1,bbeta=2,ggamma11=1,ggamma=(C[1]+aalpha)*(C[1]+bbeta)/((cx[1]+aalpha1)*(cy[1]+bbeta1))){
  ccxy=cumsum(cxy)
  ccx=cumsum(cx)
  ccy=cumsum(cy)
  cC=cumsum(C)
  EIC=log2(((ccxy+ggamma11)*(cC+aalpha)*(cC+bbeta))/((cC+ggamma)*(ccx+aalpha1)*(ccy+bbeta1)))
  VIC=((cC-ccxy+ggamma-ggamma11)/((ccxy+ggamma11)*(1+cC+ggamma)) + (cC-ccx+aalpha-aalpha1)/((ccx+aalpha1)*(1+cC+aalpha)) + (cC-ccy+bbeta-bbeta1)/((ccy+bbeta1)*(1+cC+bbeta)))/(log(2)**2)
  return(cbind(expectation=EIC,variance=VIC))
}

#examples prednisone+renal agenesia
# vC=c(5018,16000)
# vxy=c(2,3)
# vx=c(8,15)
# vy=c(101,497)
# 
# bcpnn=BCPNN(vxy,vx,vy,vC)
# ggplot()+geom_point(aes(x=c("SIDPNG","SIDPNG+1"),y=as.data.frame(bcpnn)$expectation))+
#   geom_errorbar(aes(x=c("SIDPNG","SIDPNG+1"),ymin=as.data.frame(bcpnn)$expectation-1.96*as.data.frame(bcpnn)$variance,ymax=as.data.frame(bcpnn)$expectation+1.96*as.data.frame(bcpnn)$variance))+geom_hline(aes(yintercept=0),col="red",lty=2)+
#   labs(y="IC",x="database",title="IC with BCPNN for prednisone/renal agenesia")
