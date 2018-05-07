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
DYSPLASIE3=DYSPLASIE3[-which(DYSPLASIE3=="")]

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
todel=which(apply(count.dyspatc,1,sum)<=2)
count.sym=sapply(c(DYSPLASIE3,"NONE",ATC5), function(dysp){ #le tableau de contingence
  sapply(c(DYSPLASIE3,"NONE",ATC5), function(atc) sum(grepl(dysp,data$DYSPLASIES)&grepl(atc,data$ATC)) )
})

library(igraph)
graphe=graph_from_adjacency_matrix(count.sym,weighted = T,mode = "undirected")
plot.igraph(graphe,width=E(graphe)$weight,label.cex=0.01,xlim=c(-0.01,0.01))

plot.igraph(graphe,)

compcon=decompose(graphe)
mincut=min_cut(graphe,value.only = F)


freq.atc=count.dyspatc/apply(count.dyspatc, 1, sum)
freq.dysp= t(t(count.dyspatc)/apply(count.dyspatc, 2, sum))
plot(hclust(dist(freq.atc)))
plot(hclust(dist(t(freq.dysp))))
order.atc=hclust(dist(freq.atc))$order
order.dysp=hclust(dist(t(freq.dysp)))$order
freq.order.atc=freq.atc[order.atc,order.dysp]
freq.pv=fisher.dyspatc[order.atc,order.dysp]
melt.freq.order=melt(freq.order.atc)
ggplot()+geom_raster(data=melt.freq.order,aes(x=Var2,y=Var1,fill=value))


heatmap(freq.dyspatc)

plot.phylo(as.phylo(hclust(dist(atc))),type="fan")



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


