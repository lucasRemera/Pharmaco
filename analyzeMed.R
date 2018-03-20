##########################################################
# Outils pour l'analyse de donnees Medicaments/Dysplasie #
##########################################################

#0. Chargement des donnees
setwd("~/Pharmaco") #ou on travaille -> a modifier selon le dossier de travail "working directory"
data=read.csv2("MED_DYSPLASIE.csv") # lire le fichier .csv, qui doit etre dans le meme dossier
head(data) #voir le haut de la matrice
#si erreur, esssayer 
# data=read.csv2("MED_DYSPLASIE.csv", sep=";" ) ou ... sep="," cad comment sont definis les separateur dans le fichiers csv

data$MALF3=substr(data$CODE_DYSP,1,3) #recupere les 3 premiers caracteres du code de dysplasie
data$ATC5=substr(data$CODE_ATC,1,5) #recupere les 5 premiers caracteres du code ATC

#1. Analyse univariee
#1.1 test de l'independance par malf
byMalf=sapply(sort(unique(data$MALF3)), function(j){
  fisher.test(table(data$ATC5,data$MALF3==j),simulate.p.value = T,B = 10000)$p.value})

#1.1.1 correction de bonferroni
which(byMalf<0.05/length(byMalf))

#a priori on trouve rien : faire sans correction 

#voir aussi par medicaments ...

#1.2 pour chaque combinaisons malf/med
m.pv=sapply(sort(unique(data$ATC5)), function(i){
  sapply(sort(unique(data$MALF3)), function(j){
    fisher.test(table(data$ATC5==i,data$MALF3==j),alternative = "greater")$p.value
  })
})

#2 possiblites : soit on a pour hypothese alternatives odds-ratio>1 (alternative="greater")
#soit on a H1 OR!=1 (alternative="two.sided"), il faut alors aussi recuperer l'odds-ratio et filtrer dessus

library(reshape2)
melt.pv=melt(m.pv) #restructure la matrice 

melt.pv[which(melt.pv$value<0.05),]
#1.2.2 differentes corrections
#bonferroni

alpha.bonferroni=0.05/nrow(melt.pv)

melt.pv[which(melt.pv$value<alpha.bonferroni),]

melt.pv.sort=melt.pv[order(melt.pv$value),]
#holm 
ntest=(nrow(m.pv)-1)*(ncol(m.pv)-1)
pv.holm=sapply(1:nrow(melt.pv.sort), function(k) melt.pv.sort$value[k]*(ntest+1-k))
melt.pv.sort[1:(which(pv.holm>0.05)[1]-1),]

#benjamini
pv.benjamini=sapply(1:nrow(melt.pv.sort), function(k) melt.pv.sort$value[k]*(ntest)/k)
melt.pv.sort[1:(which(pv.benjamini>0.05)[1]-1),]

##en theorie : du plus stringent au moins stringent pour les trois tests

#2. analyse multivariee / data mining
#2.1 AFC
#idee : on represente tout le tableau de contingence, dans un minimum de dimension
df=as.data.frame.matrix(table(data$MALF3,data$ATC5))
library(ade4)
afc=dudi.coa(data.frame(df))
8
scatter(afc,posieig = "none") #la proximite denote une association
atc=afc$c1
dysp=afc$l1
colnames(dysp)=paste("CS",1:8,sep = "")

library(ggplot2)
library(ggrepel)
#une representation "personnalisee" et plus souple
ggplot()+geom_text_repel(aes(x=atc$CS1,y=atc$CS2,label=rownames(atc)),col="blue")+
  geom_text_repel(aes(x=dysp$CS1,y=dysp$CS2,label=rownames(dysp)),col="red")+scale_x_continuous(limits=c(-5,5))+scale_y_continuous(limits=c(-5,5))

#on peut zoomer en modifier les x (e.g., xlimits=c(-2.5,5)) et les y (e.g., limits=c(-4,2.5))
summary(afc) #attention, ce plan ne represente qu'une partie de l'info


df.atc.dysp=rbind(atc,dysp) #concatenation des atc et dysp

library(cluster)
gap=clusGap(df.atc.dysp,kmeans,30,B=500)
nc=maxSE(gap$Tab[,3],gap$Tab[,4],method = "Tibs2001SEmax") #le nombre de clusters
print(nc)
kmeans(df.atc.dysp,nc)

#clustering hierarchique
hc=hclust(dist(df.atc.dysp))
library(ape)
ape.hc=as.phylo(hc)
plot(hc)
plot(ape.hc,cex=.8,use.edge.length = T,type="fan")

#2.2 ACM
data.acm=t(sapply(unique(data$NUM_DOSSIER), function(nodos){
  c(unique(data$MALF3),unique(data$ATC5))%in%c(data$ATC5[which(data$NUM_DOSSIER==nodos)],data$MALF3[which(data$NUM_DOSSIER==nodos)])
  #sapply(c(unique(data$MALF3),unique(data$ATC5)), function)
}))

colnames(data.acm)=c(unique(data$MALF3),unique(data$ATC5))
rownames(data.acm)=unique(data$NUM_DOSSIER)
data.acm=data.frame(apply(data.acm,2,factor))
library(ade4)
res.acm=dudi.acm(data.acm,scannf = F,nf=186)

#scatter(res.acm)
acm.value=res.acm$c1
acm.value$type=c(rep("DYSP",202),rep("MED",388-202))
acm.value$TF=sapply(strsplit(rownames(acm.value),split='[.]'),function(i) i[2])
acm.value$nom=sapply(strsplit(rownames(acm.value),split='[.]'),function(i) i[1])
library(ggplot2)
library(ggrepel)

ggplot(data=acm.value[which(acm.value$TF=="TRUE"),])+geom_point(aes(x=CS1,y=CS2,col=type))
ggplot(data=acm.value[which(acm.value$TF=="TRUE"),])+geom_text_repel(aes(x=CS1,y=CS2,label=nom,col=type))

#re-cluster
df.acm.true=acm.value[which(acm.value$TF=="TRUE"),]
hc=hclust(dist(df.acm.true[,1:186]))
library(ape)
ape.hc=as.phylo(hc)
ape.hc$tip.label=df.acm.true$nom
plot(ape.hc,tip.color = c("red","blue")[as.integer(as.factor(df.acm.true$type))],cex=.8,use.edge.length = T,type="fan")


#2.3 heatmap
data.table=table(data$ATC5,data$MALF3)

as.integer(as.factor(substr(colnames(data.table),1,2)))
colcol=substr(colnames(data.table),1,1)
colcol[which(colcol=="Q")]=substr(colnames(data.table),1,2)[which(colcol=="Q")]
colcol=as.integer(as.factor(colcol))

heatmap(t(scale(t(data.table))),RowSideColors = rainbow(13)[as.integer(as.factor(substr(rownames(data.table),1,1)))], ColSideColors = rainbow(17)[colcol])



hc.atc=(hclust(dist(t(scale(t(data.table))))))
nc=5
clus.atc=cutree(hc.atc,nc)

table(substr(rownames(data.table),1,1),clus.atc)

plot(as.phylo(hc.atc),cex=0.8,use.edge.length = T,tip.color = rainbow(nc)[clus.atc])


#work in progress

hcut<-function(x,k){
  cutree(hclust(dist(x)),k)
}

hclusCut <- function(x, k, d.meth = "euclidean", ...)
  list(cluster = cutree(hclust(dist(x, method=d.meth), ...), k=k))

clusGap(t(scale(t(data.table))),hclusCut,20,B=500)

testK=2:60
silK=sapply(testK, function(k) mean(silhouette(cutree(hc.atc,k),dist(t(scale(t(data.table)))))[,3]))
plot(testK,silK)

#library(NbClust)

#NbClust(t(scale(t(data.table))),method = "complete")
nc=4
scaled.data.table=t(scale(t(data.table)))
silhouetteMedDiscrByDysp=sapply(1:ncol(data.table), function(cc) silhouette(cutree(hc.atc,nc),dist(scaled.data.table[,cc]))[,3])
