# Script to transform tables (type: Lookup Table) from Dinamica-EGO 3.0 to R-Project and use in: Hierarchical Cluster Analysis (HCA method: Ward) defining clusters with function Ward (k-best),  optimal number of groups, saving graphs and new table
# Coding: G1S8B1TS # IGC/UFMG 2015-2016 (CC) #CreativeCommons

## Instalacion de los paquetes necesarios (Necesita elevar el nivel del usuario como /Administrador)
install.packages("dendextend", "cluster")

## Selecionar el directorio de la tabla (.csv) creada por Dinamica-EGO
setwd("C:/KDE/Endemismo/resultados/")
tabela_maxima_distancia<-read.csv("tabela_maxima_distancia.csv", sep=",", check.names=F)

## Preparar la Lookup Table para R-Project
tabela_maxima_distancia<-tabela_maxima_distancia[,c(1:2)] 
colnames(tabela_maxima_distancia)<-c("SID", "Value")
row.names(tabela_maxima_distancia)<-tabela_maxima_distancia$SID

## Matriz de distancia Ward para el análisis de Cluster
x<-tabela_maxima_distancia[,1]
xdis<-dist(tabela_maxima_distancia, method="euclidean")

## Cluster usando Ward
dist.euc.ward <- hclust(xdis, method="ward.D")
plot(dist.euc.ward,main="Dist - Euclidean - Ward")

library(cluster)
## Guarda el grafico numero optimo de cluster
jpeg("Silhouette-optimal number of clusters.jpeg", res=300, width=2000, height=1000)
asw <- numeric(length(x))
for (k in 2:(length(x)-1)) {
  sil <- silhouette(cutree(dist.euc.ward, k=k), xdis)
  asw[k] <- summary(sil)$avg.width
}
k.best <- which.max(asw)
op<-par(no.readonly=TRUE)
plot(1:length(x), asw, type="h", xlab="k (number of groups)", ylab="Average width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="red", font=2,
     col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
dev.off()

## Guarda el grafico dendrograma
library(dendextend)
jpeg("Dendrogram Euclidean distance, Ward method.jpeg", res=300, width=2000, height=1000)
dend<-as.dendrogram(dist.euc.ward)
dend <- color_branches(dend, k=k.best, col=topo.colors)
dend <- hang.dendrogram(dend,hang_height=0.2)
dend <- set(dend, "labels_cex", 0.5)
par(mar = c(3,3,3,3))
plot(dend, horiz = F, nodePar = list(cex = .07))
dev.off()

##  Crea y guarda una tabla de grupos basado en el K(Numero optimo)
groups<-as.data.frame(sort(cutree(dend, k=k.best)))
groups$Key<-as.numeric(row.names(groups))
row.names(groups)<-NULL
groups$Group<-groups[,1]
groups<-groups[,c(2,3)]
colnames(groups)<-c("Key","Group")
write.csv(groups, "tabela_grupos.csv",row.names = F)

