BiocManager::install('fibroEset')
library(fibroEset)
library(stats)
library(MASS)
data(fibroEset)
ann <- fibroEset$species
dat <- exprs(fibroEset)

dat <-dat[sample(nrow(dat),50),]
colnames(dat) <- paste(colnames(dat),ann,' ')

hc<-hclust(dist(t(dat),method= "manhattan"), method = "median")

plot(hc,
     main = 'Hierarchical Clustering Dendrogram among 50 genes \n firbroEset HGU95Av2 data',
     xlab = 'Sample\n h: human (Homo sapiens), b: bonobo (Pan paniscus), g: gorilla (Gorilla gorilla)'
     ,sub = '')

hc
hm.rg <-   c("#FF0000","#CC0000","#990000","#660000","#330000","#000000","#000000","#0A3300","#146600","#1F9900","#29CC00","#33FF00")
heatmap(dat,main ='Heatmap among 50 genes in firbroEset HGU95Av2 data', col=hm.rg,
        xlab = 'Sample\n h: human (Homo sapiens), b: bonobo (Pan paniscus), g: gorilla (Gorilla gorilla)',
        ylab = 'probset')

pca.dat<- prcomp(t(dat),center = T,scale. = T,cor=F)$x[,1:2]
kc <- kmeans(pca.dat, centers=3,iter.max = 20)
plot(pca.dat, col = kc$cluster,cex=1,main='K-means Clustering PCA Scatter Plot \n k=3',lwd=1.5)
points(kc$centers,col=1:3,pch='*',cex=2.5)

