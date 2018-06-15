setwd("~/Pharmaco") #ou on travaille -> a modifier selon le dossier de travail "working directory"
#data=read.csv2("DYSP_ATC.csv") # lire le fichier .csv, qui doit etre dans le meme dossier
data=read.csv2("DYSP_ATC1.csv")

# head(data)
# DYSPLASIES                                                             ATC
# 1       K07.0|K07.0|K07.1|K07.1|Q37.8|Q37.8|Q87.08|Q87.08 N02AE01|N05BA12|N02AE01|N05BA12|N02AE01|N05BA12|N02AE01|N05BA12
# 2 P29.2|Q05.9|Q06.1|Q21.0|Q21.1|Q21.11|Q23.10|Q23.4|Q25.4                                                        ||||||||
# 3                                             Q05.7|Q06.1                                                               |
# 4                                 P29.2|P94.2|Q21.1|Q90.0                                                             |||

#une entree pour chaque cas, une colonne avec tous les codes de dysplasies associées à ce cas et une autre avec tous les codes ATC associés à ce cas

#vecteur des 3 premiers caractères du code dysplasie
DYSPLASIE=unique(unlist(strsplit(as.character(data$DYSPLASIES),split = "[|]")))
DYSPLASIE3=unique(substr(DYSPLASIE,1,3))
#DYSPLASIE3=DYSPLASIE3[-which(DYSPLASIE3=="")]

#vecteur des 5 premiers caractères du code ATC
ATC=unique(unlist(strsplit(as.character(data$ATC),split = "[|]")))
ATC=ATC[-which(ATC=="")]
ATC5=unique(substr(ATC,1,5))

#conversion
data=data.frame(apply(data, 2, as.character),stringsAsFactors = FALSE)

#creation du code "NONE" quand il n'y a pas d'ATC
data$ATC[which(apply(sapply(LETTERS,function(lett) grepl(lett,data$ATC)),1,sum)==0)]="NONE"

# table(grepl(DYSPLASIE3[1],data$DYSPLASIES),grepl(ATC5[1],data$ATC))
# 
# table(grepl("Q86",data$DYSPLASIES),grepl("C09CA",data$ATC))

##1.Analyse descriptive multivariee : AFC (cf. wiki)
count.dyspatc=sapply(DYSPLASIE3, function(dysp){ #le tableau de contingence
  sapply(c("NONE",ATC5), function(atc) sum(grepl(dysp,data$DYSPLASIES)&grepl(atc,data$ATC)) )
})

library(ade4)
afc=dudi.coa(data.frame((count.dyspatc)))
11
summary(afc)
#colonne Cumulative Projected inertia, Ax1:2 pour connaitre l'information contenue dans le premier plan
scatter(afc,posieig = "none")
atc=afc$l1
dysp=afc$c1
colnames(dysp)=paste("CS",1:11,sep = "")
colnames(atc)=paste("CS",1:11,sep = "")

#1.2 représentation des codes ATC
ggplot()+geom_text_repel(aes(x=atc$CS1,y=atc$CS2,label=rownames(atc)),col=c("red",rep("blue",nrow(atc)-1)))

#1.3 ATC+Dysplasies
#xlim et ylim pour zoomer
ggplot()+geom_text_repel(aes(x=atc$CS1,y=atc$CS2,label=rownames(atc)),col=c("red",rep("blue",nrow(atc)-1)))+
  geom_text_repel(aes(x=dysp$CS1,y=dysp$CS2,label=rownames(dysp)),col="green")+
  xlim(-5,7.5)+ylim(-2,4)

#1.4 : quels sont les codes ATC avec une repartition parmi les dysplasies la plus anormale?
distanceANONE=apply(atc,1,function(i) sqrt(sum((i-atc[1,])**2) ))
sort(distanceANONE,decreasing = T)

#1.5 : quels sont les dysplasies avec une repartition parmi les ATC la plus anormale?
distanceAZero=apply(dysp,1,function(i) sqrt(sum((i)**2) ))
sort(distanceAZero,decreasing = T)


##2. tests de Fisher paire par paire
fisher.dyspatc=sapply(DYSPLASIE3, function(dysp){
  sapply(ATC5, function(atc) fisher.test(table(grepl(dysp,data$DYSPLASIES),grepl(atc,data$ATC)),alternative = "greater")$p.value  )
})

melt.dyspatc=melt(fisher.dyspatc)
melt.dyspatc[Benjamini(melt.dyspatc$value)$validate,] #selection apres correction de Benjamini
melt.dyspatc[which(melt.dyspatc$value<Benjamini(melt.dyspatc$value[-which(melt.dyspatc$value==1)])$alphaCorrige),]

##alternative pour etre moins stringent
# todel=which(apply(count.dyspatc,1,sum)<=2)
# count.sym=sapply(c(DYSPLASIE3,"NONE",ATC5), function(dysp){ #le tableau de contingence
#   sapply(c(DYSPLASIE3,"NONE",ATC5), function(atc) sum(grepl(dysp,data$DYSPLASIES)&grepl(atc,data$ATC)) )
# })

#count.dyspatc=count.dyspatc[-1,]
# natc=nrow(count.dyspatc)
# ndysp=ncol(count.dyspatc)
# count.sym2=matrix(0,ncol=sum(dim(count.dyspatc)),nrow=sum(dim(count.dyspatc)))
# 
# count.sym2[(ndysp+1):(ndysp+natc),1:ndysp]=count.dyspatc
# nom=c(colnames(count.dyspatc),rownames(count.dyspatc))
# colnames(count.sym2)=nom
# rownames(count.sym2)=nom
# 
# #count.sym2=count.sym2+t(count.sym2)
# (spdysp=spectralPartitionnement(count.sym2,norme = F,k=2))
# nom[which(spdysp==1)]
# 
# library(igraph)
# graphe=graph_from_adjacency_matrix(count.sym,weighted = T,mode = "undirected")
# plot.igraph(graphe,width=E(graphe)$weight,label.cex=0.01,xlim=c(-0.01,0.01))
# 
# plot.igraph(graphe)
# 
# compcon=decompose(graphe)
# mincut=min_cut(graphe,value.only = F)
# 
# 
# freq.atc=count.dyspatc/apply(count.dyspatc, 1, sum)
# freq.dysp= t(t(count.dyspatc)/apply(count.dyspatc, 2, sum))
# plot(hclust(dist(freq.atc)))
# plot(hclust(dist(t(freq.dysp))))
# order.atc=hclust(dist(freq.atc))$order
# order.dysp=hclust(dist(t(freq.dysp)))$order
# freq.order.atc=freq.atc[order.atc,order.dysp]
# freq.pv=fisher.dyspatc[order.atc,order.dysp]
# melt.freq.order=melt(freq.order.atc)
# ggplot()+geom_raster(data=melt.freq.order,aes(x=Var2,y=Var1,fill=value))
# 
# 
# heatmap(freq.dyspatc)
# 
# plot.phylo(as.phylo(hclust(dist(atc))),type="fan")
# 
# 
# 
# df.atc.dysp=rbind(atc,dysp) #concatenation des atc et dysp
# 
# library(cluster)
# gap=clusGap(df.atc.dysp,kmeans,30,B=500)
# nc=maxSE(gap$Tab[,3],gap$Tab[,4],method = "Tibs2001SEmax") #le nombre de clusters
# print(nc)
# kmeans(df.atc.dysp,nc)
# 
# #clustering hierarchique
# hc=hclust(dist(df.atc.dysp))
# library(ape)
# ape.hc=as.phylo(hc)
# plot(hc)
# plot(ape.hc,cex=.8,use.edge.length = T,type="fan")

## ACM ##
ATC_=ATC5[-which(nchar(ATC5)<5)] 
TABLE_DYSP=sapply(DYSPLASIE3, function(dysp) as.factor(grepl(dysp,data$DYSPLASIES)))
TABLE_ATC=sapply(c("NONE",ATC_), function(atc) as.factor(grepl(atc,data$ATC)))

burt=as.data.frame(cbind(TABLE_DYSP,TABLE_ATC))
#df.burt=as.data.frame(apply(burt,2,as.factor))
#subburt=burt[1:5000,]
#torem=which(apply(subburt,2,function(i) sum(grepl("TRUE",i)))==0)
#subburt=subburt[,-torem]

#acm=dudi.acm(subburt)
#s.label(acm$c1)
#s.label(acm$l1)

burt.acp=apply(burt,2,function(i) as.integer(as.logical(as.character((i)))))
medacp=dudi.pca(burt.acp)
s.corcircle(medacp$c1)
s.label(medacp$l1)
ggplot()+geom_point(data=medacp$l1,aes(x=RS1,y=RS2))+#xlim(0,.55)+ylim(-1,10)+
  xlim(0,0.05)+ylim(-1,10)+
  geom_text_repel(aes(x=medacp$c1$CS1,y=medacp$c1$CS2,label=rownames(medacp$c1),col=factor(nchar(rownames(medacp$c1)))))+labs(col="")


#data[which(medacp$l1[,1]>.2&medacp$l1[,1]<1),]

#fisher.test(table(data$ATC=="NONE",data$DYSPLASIES),simulate.p.value = T,B = 100000)

#######################################
# HEATMAP of standardized frequencies #
#######################################
contingence=count.dyspatc[which(rownames(count.dyspatc)=="NONE"|nchar(rownames(count.dyspatc))==5),]
# N=sum(contingence)
# meltedcontingence=expand.grid(atc=apply(contingence, 1, sum),dysp=apply(contingence, 2, sum))
# #M0=matrix(meltedcontingence[,2],nrow=nrow(contingence))
# M0=matrix(apply(meltedcontingence, 1, prod)/(N),nrow=nrow(contingence))
# fij=matrix(apply(meltedcontingence, 1, prod)/(N**2),nrow=nrow(contingence))  
# zij=(contingence-M0)/(1.96*sqrt(fij*(1-fij)*N)) #standardization
# melt.zij=melt(zij) 
# ggplot()+geom_raster(data=melt.zij,aes(x=Var1,y=Var2,fill=value))+
#   scale_fill_gradientn(colours = c("green","white","red","darkred"),values = c(0,0,1))+
#   labs(fill="Z value")+ylab("dysplasie")+xlab("ATC")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y = element_text(size=5))
# heatmap(zij)
# 
# hatc=(hclust(dist(zij)))
# hdysp=hclust(dist(t(zij)))
# 
# #ordered with hclust
# 
# contingence.order=contingence[hatc$order,hdysp$order]
# meltedcontingenceorder=expand.grid(atc=apply(contingence.order, 1, sum),dysp=apply(contingence.order, 2, sum))
# #M0=matrix(meltedcontingence[,2],nrow=nrow(contingence))
# M0order=matrix(apply(meltedcontingenceorder, 1, prod)/(N),nrow=nrow(contingence))
# fijorder=matrix(apply(meltedcontingenceorder, 1, prod)/(N**2),nrow=nrow(contingence))  
# zijorder=(contingence.order-M0order)/(1.96*sqrt(fijorder*(1-fijorder)*N)) #standardization
# melt.zijorder=melt(zijorder) 
# ggplot()+geom_raster(data=melt.zijorder,aes(x=Var1,y=Var2,fill=value))+
#   scale_fill_gradientn(colours = c("green","white","red","darkred"),values = c(0,0,1))+
#   labs(fill="Z value")+ylab("dysplasie")+xlab("ATC")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y = element_text(size=5))
# 
# ggplot()+geom_raster(data=melt.zijorder,aes(x=Var1,y=Var2,fill=ifelse(value>0,value,NA)))+
#   scale_fill_gradient2(low="blue",mid="pink",high="darkred",midpoint=1,na.value = "white")+
#   labs(fill="Z value")+ylab("dysplasie")+xlab("ATC")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y = element_text(size=5))

standardizeMatrix=function(m,replaceNA=TRUE,phialpha=1.96){
  N=sum(m)
  meltedcontingence=expand.grid(atc=apply(m, 1, sum),dysp=apply(m, 2, sum))
  #M0=matrix(meltedcontingence[,2],nrow=nrow(contingence))
  M0=matrix(apply(meltedcontingence, 1, prod)/(N),nrow=nrow(m))
  fij=matrix(apply(meltedcontingence, 1, prod)/(N**2),nrow=nrow(m))  
  zij=(m-M0)/(phialpha*sqrt(fij*(1-fij)*N)) #standardization
  if(replaceNA) zij[is.na(zij)]<-0
  return(zij)
}

plotHeatmap=function(zij,toOrder=TRUE){
  nr=nrow(zij)
  nc=ncol(zij)
  if(toOrder){
    hatc=(hclust(dist(zij)))
    hdysp=hclust(dist(t(zij)))
    zij=zij[hatc$order,hdysp$order]
  }
  melt.zij=melt(zij) 
  ggplot()+geom_raster(data=melt.zij,aes(x=Var1,y=Var2,fill=value))+
    scale_fill_gradientn(colours = c("green","white","red","darkred"),values = c(0,0,1))+
    labs(fill="Z value")+ylab("dysplasie")+xlab("ATC")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5),axis.text.y = element_text(size=5))
  
}

plotHeatmap(standardizeMatrix(contingence))

###############
# simulations #
###############

associationH0=function(fij,N){
  apply(fij,1,function(i) rbinom(ncol(fij),N,i))  
}

pairFisher=function(contingence){
  nc=ncol(contingence)
  nr=nrow(contingence)
  mf=sapply(1:nr, function(rr){
    sapply(1:nc, function(cc){
      em=contingence[rr,cc]
      nem=sum(contingence[-rr,cc])
      enm=sum(contingence[rr,-cc])
      nenm=sum(contingence[-rr,-cc])
      if(em==0|em+enm==0|nem+em==0) return(1)
      else return(fisher.test(matrix(c(em,enm,nem,nenm),nrow=2),alternative = "greater" )$p.value )
    })
  })
  mf=t(mf)
  colnames(mf)=colnames(contingence)
  rownames(mf)=rownames(contingence)
  return((mf))
}


simulateSpecifity=function(fij,Nobs,Nsim=100){
  return(sapply(1:Nsim, function(ii){
    m0=associationH0(fij,Nobs)
    return(melt(pairFisher(m0))[,3])
  }))
}

specificity=function(pvalues,alphas=seq(0,1,by=0.01)){
  x=sapply(alphas,function(a) sum(pvalues<=a))/length(pvalues)
  names(x)=alphas
  return(x)
}


z0=standardizeMatrix(associationH0(fij,N))
plotHeatmap(z0)

dummym=matrix(rpois(6,10),nrow=2)
pairFisher(dummym)
fishC=pairFisher(contingence)

dummy0=simulateSpecifity(fij,N,100)
#plot(sort(dummy0[,1]))
alphaseuil=seq(0,0.1,by=0.001)
spe0=apply(dummy0,2,function(ii) specificity(ii,alphaseuil))
spe1=specificity(as.numeric(fishC),alphaseuil)

ggplot()+geom_line(data=melt(spe0),aes(x=Var1,y=value,group=Var2,col=factor(Var2)),lty=2)+
  geom_point(aes(x=alphaseuil,y=spe1),col="red")+ylim(0,0.012)+xlim(0,.05)
# cspe0=apply(spe0,2,cumsum)
# ggplot()+geom_line(data=melt(cspe0),aes(x=Var1,y=value,group=Var2,col=factor(Var2)),lty=2)+
#   geom_point(aes(x=alphaseuil,y=cumsum(spe1)),col="red")+ylim(0,0.2)+xlim(0,.05)


##############
# sensitvity #
##############
associationH1=function(fij,N,i=NULL,j=NULL,RR=2){
  rr=rep(1,ncol(fij))
  rr[j]=RR
  m=sapply(1:nrow(fij),function(ii){
    if(ii!=i) rbinom(ncol(fij),N,fij[ii,])
    else rbinom(ncol(fij),N,fij[ii,]*rr)
  } )  
  return(t(m))
}

simulateSensitivity=function(fij,Nobs,Nsim=100,i=NULL,j=NULL,RR=c(1,2),pdsRandom=TRUE){
  if(is.null(i)|is.null(j)){
    if(pdsRandom){
      wi=apply(fij, 1, sum)
      wj=apply(fij, 2, sum)
    }else{
      wi=rep(1,nrow(fij))
      wj=rep(2,ncol(fij))
    }
    wi=wi/sum(wi)
    wj=wj/sum(wj)
  }
  if(length(RR)==1) RR=rep(RR,2)

  pv=sapply(1:Nsim, function(k){
    if(is.null(i)) i=sample(1:nrow(fij),1,prob = wi)
    if(is.null(j)) j=sample(1:ncol(fij),1,prob = wj)
    mas=associationH1(fij,Nobs,i,j,runif(1,RR[1],RR[2]))
    em=mas[i,j]
    nem=sum(mas[-i,j])
    enm=sum(mas[i,-j])
    nenm=sum(mas[-i,-j])
    if(em==0|em+enm==0|nem+em==0) return(1)
    else return(fisher.test(matrix(c(em,enm,nem,nenm),nrow=2),alternative = "greater" )$p.value )
  })
  return(pv)
}

pvH1nor=simulateSensitivity(fij,N,RR=2,Nsim = 10000,pdsRandom = F)
speH1nor=specificity(pvH1nor,alphaseuil)
ggplot()+geom_line(data=melt(spe0),aes(x=Var1,y=value,group=Var2,col=factor(Var2)),lty=2)+
  geom_point(aes(x=alphaseuil,y=speH1nor),col="red")+#ylim(0,0.012)+xlim(0,.05)+
  theme(legend.position = "none")



alphaseuil=seq(0,1,by=0.001)
spe0=apply(dummy0,2,function(ii) specificity(ii,alphaseuil))
speH1nor=specificity(pvH1nor,alphaseuil)
TFP=apply(spe0,1,mean)
TVP=speH1nor
plot(TFP,TVP,type="l",xlim=c(0,.15),ylim=c(0,.15))
abline(a = 0,b=1)
plot(alphaseuil,TVP-TFP,ylim=c(-.05,.05))

pvH1r=sapply(seq(1,3,by=.1), function(rr) simulateSensitivity(fij,N,RR=rr,Nsim = 5000,pdsRandom = F))
pvH1s=apply(pvH1r,2, function(ii)  specificity(ii,alphaseuil))
mpvh1s=melt(pvH1s)

ggplot()+geom_line(data=melt(spe0),aes(x=Var1,y=value,group=Var2,linetype="TFP"))+
  geom_line(data=mpvh1s,aes(x=Var1,y=value,group=Var2,col=(Var2/20+1),linetype="TVP"))+ylim(0,0.06)+xlim(0,.08)+
  labs(col="RR",linetype="")+
  scale_linetype_manual(c("TFP","TVP"),values=c(2,1))+
  xlab("p-value threshold")+
  ylab("positive rate")+
  ggtitle("True and false positive rate")
  #theme(legend.position = "none")


rTP=apply(pvH1s, 2, function(i) i/TFP)
ggplot()+geom_line(data=melt(rTP),aes(x=Var1,y=value,group=Var2,col=(Var2-1)/10+1))+labs(col="RR")+
  ylim(0.9,10)+xlim(0,.05)

###########################################
# influence de la taille de l'echantillon #
###########################################
vn=seq(5000,5000,by=10000)
dummy0N=sapply(vn, function(nn) simulateSpecifity(fij,nn,100)) 
