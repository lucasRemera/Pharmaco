sidp1=read.csv2("extraction sidp1.csv",header = T)
sidpng=read.csv2("extraction confusion F et M.csv",header = T)
annee1=as.numeric(as.character(ifelse(nchar(sidp1$NODOS)==6,substr(sidp1$NODOS,1,1),substr(sidp1$NODOS,1,2))))
anneeng=as.numeric(as.character(substr(sidpng$NUM_DOSSIER,3,4)))
q601=grepl("Q60",sidp1$DYSPLASIES)
q60ng=grepl("Q60",sidpng$DYSPLASIES)
h02ab1=grepl("H02AB07",sidp1$ATC)
h02abng=grepl("MC",sidpng$ATC) 
annee=c(annee1,anneeng)
q60=as.numeric(c(q601,q60ng))
prednisone=as.numeric(c(h02ab1,h02abng))
uannee=unique(annee)
vC=as.numeric(table(annee))
vx=as.integer(table(prednisone,annee)[2,])
vy=as.integer(table(q60,annee)[2,])
vxy=table(prednisone,q60,annee)[2,2,]
bcpnna=BCPNN(vxy,vx,vy,vC,aalpha = 1,bbeta = 1,ggamma = (length(q60)**2)/(sum(vx)*sum(vy)))
bcpnna=BCPNN(vxy,vx,vy,vC,aalpha = 1,bbeta = 1)
bcpnna=BCPNN(vxy,vx,vy,vC)
ggplot()+geom_point(aes(x=as.factor(2000+uannee),y=as.data.frame(bcpnna)$expectation))+
     geom_errorbar(aes(x=as.factor(2000+uannee),ymin=as.data.frame(bcpnna)$expectation-1.96*as.data.frame(bcpnna)$variance,ymax=as.data.frame(bcpnna)$expectation+1.96*as.data.frame(bcpnna)$variance))+geom_hline(aes(yintercept=0),col="red",lty=2)+
     labs(y="IC",x="year",title="IC with BCPNN for prednisone/renal agenesia")

qm=data.frame(cbind(a=2000+uannee,e=vx-vxy,em=vxy,m=vy-vxy))
qmm=melt(qm,id.vars = c(1))
ggplot()+geom_bar(aes(factor(qmm$a),qmm$value,fill=qmm$variable),stat="identity")


#######################
# avec nouveaux codes #
#######################
code=read.csv2("nouveaux_codes.csv")
code1=sapply(as.character(sidp1$ATC), function(vatc){
  #print(vatc)
  if(sum(sapply(LETTERS, function(l) grepl(l,vatc)))==0) return("none")
  else{
    uatc=unique(strsplit(vatc,split = "[|]")[[1]])
    #return(which(code$ATC %in% uatc))
    return(paste(unique(code$CODES[which(code$ATC %in% uatc)]),collapse = "|"))
  }
})


codeng=sapply(as.character(sidpng$ATC), function(vatc){
  #print(vatc)
  if(sum(sapply(LETTERS, function(l) grepl(l,vatc)))==0) return("none")
  else{
    return(paste(unique(strsplit(vatc,split = "[|]")[[1]]),collapse = "|"))
    
  }
})

codes=c(code1,codeng)
ucode=c("none",as.character(na.omit(unique(code$CODES))))
#ucode=ucode[-which(ucode=="")]

mcode=sapply(ucode, function(a) as.numeric(grepl(a,codes)))

nnn=apply(mcode,2,sum)
ucoden=ucode[which(nnn>0)]
mcoden=mcode[,which(nnn>0)]
lasso=cv.glmnet(mcoden,as.factor(q60),family = "binomial",alpha=1,lambda = 10**seq(-6,4,by=.2))
lasso$glmnet.fit$beta[,which(lasso$lambda== lasso$lambda.min)]

lasso0=cv.glmnet(mcoden,as.factor(q60),family = "binomial",alpha=1,lambda = 10**seq(-6,4,by=.2))

plot(lasso)
lasso0$glmnet.fit$beta[,which(lasso0$lambda==lasso0$lambda.min)]

#
mmed=mcoden[,names(which(lasso$glmnet.fit$beta[,33]>0))]




imc=ifelse(as.integer(as.character(sidpng$IMC))>0,as.integer(as.character(sidpng$IMC)),NA)
hta=ifelse(sidpng$HYPERTENSION=="OUI","OUI","NON")
rang=as.integer(as.character(sidpng$NB_VIVANTS))
terme=as.integer(as.character(sidpng$TERME_SEMAINES))
tabac=ifelse(sidpng$TABAC%in%c("0","9","1","NULL"),"0",as.character(as.integer(as.character(sidpng$TABAC))-1))
diab=ifelse(sidpng$DIABETE=="OUI","OUI","NON")
alc=ifelse(sidpng$ALCOOL%in%c("0","9","1","NULL"),"0",as.character(as.integer(as.character(sidpng$ALCOOL))-1))
# 
# #alc=ifelse(!Donneesconfusion2$ALCOOL%in%c("2","3"),"NON","OUI")
drog=ifelse(sidpng$DROGUE=="OUI","OUI","NON")
age=ifelse(as.integer(as.character(sidpng$AGE_NAISSANCE))>0,as.integer(as.character(sidpng$AGE_NAISSANCE)),NA)
inf=ifelse(sidpng$INFECTION_UN_TRI=="1","OUI","NON")
# ml2=glm(z~med2+age+diab+imc+tabac+alc+drog,family = "binomial")
# summary(ml2)
# 
# #sidp1
# sidp1=read.csv2("extraction sidp1.csv",header = T)
#med1=grepl("H02AB07",sidp1$ATC)
# malf1=as.integer(grepl("Q60",sidp1$DYSPLASIES))
age1=ifelse(sidp1$AGEMER!="NULL",as.integer(as.character(sidp1$AGEMER)),NA)
diab1=ifelse(sidp1$DIABETE=="O","OUI","NON")
imc1=as.integer(as.character(sidp1$POIMER))/((0.01*as.integer(as.character(sidp1$TAILMER)))**2)
tabac1=ifelse(as.character(sidp1$MERTABA)=="NULL","0",as.character(sidp1$MERTABA))
alc1=ifelse(as.character(sidp1$MERALCO)=="NULL","0",as.character(sidp1$MERALCO))
drog1=ifelse(sidp1$CONSODROG=="O","OUI","NON")
# sidp1$MERINF1T
inf1=ifelse(sidp1$MERINF1T=="1","OUI","NON")

AGE=c(age1,age)
DIAB=c(diab1,diab)
TABAC=c(tabac1,tabac)
attach(data.frame(mmed))

#y=cbind(mmed,AGE,DIAB,TABAC)
#y[,8]=as.numeric(y[,8])
glm1=glm(q60~AGE+DIAB+TABAC+HR+JW+FQ+PF+MC+LH+IU,family="binomial")


dfd=data.frame(AGE,DIAB,TABAC,HR,JW,FQ,PF,MC,LH,IU)


xfactors <- model.matrix(q60 ~ DIAB+TABAC+HR+JW+FQ+PF+MC+LH+IU)[, -1]
x        <- as.matrix(data.frame(AGE,  xfactors))
lasso2=cv.glmnet(x,as.factor(q60),family = "binomial",alpha=1,lambda = 10**seq(-6,4,by=.2))
