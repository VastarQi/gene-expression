BiocManager::install('edgeR')
library(corrplot)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(pcaMethods)
library(edgeR)
library(fibroEset)
library(stats)
library(MASS)

dat <- read.table("C:/Users/QiHY///OneDrive/JHU/gene expression/final/GSE38941_series_matrix.txt",row.names = 1, header=T,comment.char = '!' )
dat <- 2^dat
dat.nor <- dat[,1:10]
dat.ALF <- dat[,11:27]


pearsonmatrix <- cor(dat,method = 'pearson',use = 'pairwise.complete.obs')
col <- colorRampPalette(c('blue','white','red'))(20)
p <- pheatmap(pearsonmatrix,col=col,clustering_distance_rows='correlation',clustering_distance_cols='correlation',border=F,
              main = 'Correlation heat map among the arrays')
p

dat.mean <- apply(log2(dat),2,mean)
dat.sd <- sqrt(apply(log2(dat),2,var)) 
dat.cv <- dat.sd/dat.mean 
plot(dat.mean,dat.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21,cex=1.1)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.7)


dat.avg <- apply(pearsonmatrix,1,mean)
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avgr",main="Average correlation plot",cex.main=1.5,axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.2)
text(dat.avg,label=dimnames(dat)[[2]],pos=1,cex=0.7)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

dat.gene.mean <- apply(log2(dat),1,mean)
dat.gene.sd <- sqrt(apply(log2(dat),1,var)) 
dat.gene.cv <- dat.gene.sd/dat.gene.mean 
plot(dat.gene.mean,dat.gene.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.gene.mean,dat.gene.cv,bg="lightblue",col=1,pch=21,cex=0.5)


#visualization of outlier nor

pearsonmatrix <- cor(dat.nor,method = 'pearson',use = 'pairwise.complete.obs')
col <- colorRampPalette(c('blue','white','red'))(20)
p <- pheatmap(pearsonmatrix,col=col,clustering_distance_rows='correlation',clustering_distance_cols='correlation',border=F,
              main = 'Correlation heat map among the arrays')
p

dat.mean <- apply(log2(dat.nor),2,mean)
dat.sd <- sqrt(apply(log2(dat.nor),2,var)) 
dat.cv <- dat.sd/dat.mean 
plot(dat.mean,dat.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21,cex=1.1)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.7)


dat_scaled <- scale(t(dat.nor)[,-1])
hc<-hclust(dist(dat_scaled,method = "euclidean"),method = "ward.D2")
hcd <- plot(hc,hang = -0.01,cex=0.7,main='')


dat.avg <- apply(pearsonmatrix,1,mean)
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avgr",main="Average correlation plot",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.2)
text(dat.avg,label=dimnames(dat.nor)[[2]],pos=1,cex=0.7)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat.nor)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


dat.gene.mean <- apply(log2(dat.nor),1,mean)
dat.gene.sd <- sqrt(apply(log2(dat.nor),1,var)) 
dat.gene.cv <- dat.gene.sd/dat.gene.mean 
plot(dat.gene.mean,dat.gene.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.gene.mean,dat.gene.cv,bg="lightblue",col=1,pch=21,cex=0.5)

# ALF
pearsonmatrix <- cor(dat.ALF,method = 'pearson',use = 'pairwise.complete.obs')
col <- colorRampPalette(c('blue','white','red'))(20)
p <- pheatmap(pearsonmatrix,col=col,clustering_distance_rows='correlation',clustering_distance_cols='correlation',border=F,
              main = 'Correlation heat map among the arrays')
p

dat.mean <- apply(log2(dat.ALF),2,mean)
dat.sd <- sqrt(apply(log2(dat.ALF),2,var)) 
dat.cv <- dat.sd/dat.mean 
plot(dat.mean,dat.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21,cex=1.1)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.7)


dat_scaled <- scale(t(dat.ALF)[,-1])
hc<-hclust(dist(dat_scaled,method = "euclidean"),method = "ward.D2")
hcd <- plot(hc,hang = -0.01,cex=0.7,main='')


dat.avg <- apply(pearsonmatrix,1,mean)
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avgr",main="Average correlation plot",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.2)
text(dat.avg,label=dimnames(dat.ALF)[[2]],pos=1,cex=0.7)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat.ALF)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")


dat.gene.mean <- apply(log2(dat.ALF),1,mean)
dat.gene.sd <- sqrt(apply(log2(dat.ALF),1,var)) 
dat.gene.cv <- dat.gene.sd/dat.gene.mean 
plot(dat.gene.mean,dat.gene.cv,main="CV vs mean plot",xlab="Mean",ylab="CV",col='blue',cex=1.5,type="n")
points(dat.gene.mean,dat.gene.cv,bg="lightblue",col=1,pch=21,cex=0.5)

# muti test

aov.all.genes <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep('A',length(x1)), rep('B',length(x2)), rep('C',length(x3)),rep('D',length(x4)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c('factor','express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)}
  
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

control.col <- c(13,15,16,17,18,19,24,29,30)
AD.col <- c(1:12,14,20:23,25:28,31)
Incipient.col <- c(2,3,13,20,26,31)
Moderate.col <- c(7:9,12,21,22,25,28)
Severe.col <- c(1,4,5,6,10,11,27)


control <- dat[,control.col]
AD <- dat[,AD.col]

t.pv <- apply(dat,1,t.test.all.genes,s1=c(1:10),s2=c(11:27))
pv
hist(t.pv,xlab ='p value of t-test' )


dat.p.adj <- p.adjust(p = t.pv,method = 'holm')
dat.p.adj
hist(dat.p.adj,xlab='adjusted p value')
length(dat.p.adj[dat.p.adj<0.01])

name <- as.vector(names(dat.p.adj[dat.p.adj<0.01]))
dat.selected <- dat[name,]



dat.selected.p <- dat.p.adj[dat.p.adj<0.01]
hist(dat.selected.p,main='Histogram of selected genes p-value distribution',xlab = 'p-value')

aov.all.genes <- function(x,s1,s2,s3,s4) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  x4 <- as.numeric(x[s4])
  fac <- c(rep('A',length(x1)), rep('B',length(x2)), rep('C',length(x3)),rep('D',length(x4)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3,x4))
  names(a.dat) <- c('factor','express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}

aov.run <- apply(dat,1,aov.all.genes,s1=control.col,s2=Incipient.col,s3=Moderate.col,s4=Severe.col)
hist(aov.run)
length(aov.run[aov.run<0.001])
aov.run.adj <- p.adjust(aov.run,method = 'BH')
hist(aov.run.adj)
length(aov.run.adj[aov.run.adj<0.01])


aov.all.genes <- function(x,s1,s2,s3) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  x3 <- as.numeric(x[s3])
  fac <- c(rep('a',length(x1)), rep('b',length(x2)), rep('c',length(x3)))
  a.dat <- data.frame(as.factor(fac),c(x1,x2,x3))
  names(a.dat) <- c('factor','express')
  p.out <- summary(aov(express~factor, a.dat))[[1]][1,5]
  #p.out <- summary(aov(express~factor, a.dat))[[1]][1,4]	# use to get F-statistic
  return(p.out)
}
aov.run <- apply(dat,1, aov.all.genes,s1=control.col,s2=Incipient.col,s3=Moderate.col)
hist(aov.run)
adj<- TukeyHSD(aov.run,conf.level=.95)

#clustering and PCA
gene.pca <- prcomp(t(dat.selected),cor=F)
dat.selected.pca.var <- round(gene.pca$sdev ^2 / sum(gene.pca$sdev^2)*100,2 )
plot(c(1:length(dat.selected.pca.var)),dat.selected.pca.var,type="b",
     xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1)
title("Scree plot showing % variability explained by each eigenvalue\nGSE38941")

gene.pca.loading <- gene.pca$x[,1:2]
plot(range(gene.pca.loading[,1]),range(gene.pca.loading[,2]),type="n",xlab='p1',ylab='p2',
     main='PCA plot of GSE38941\np2 vs. p1')
points(gene.pca.loading[,1][as.vector(colnames(dat[1:10]))], gene.pca.loading[,2][as.vector(colnames(dat[1:10]))],
       col=1,bg='red',pch=21,cex=1)
points(gene.pca.loading[,1][as.vector(colnames(dat[11:27]))], gene.pca.loading[,2][as.vector(colnames(dat[11:27]))],
       col=1,bg='blue',pch=21,cex=1)
legend('topleft',title='TYPE',c('Normal','ALF'),pch = 21,col=1,pt.bg  = c('red','blue'))

#classification

  # from selected gene
pca.dat<- prcomp(t(dat.selected),center = T,scale. = T,cor=F)$x[,1:2]
kc <- kmeans(pca.dat, centers=2,iter.max = 20)
plot(pca.dat, col = kc$cluster,cex=1,main='K-means Clustering PCA Scatter Plot \n k=2',lwd=1.5)
points(kc$centers,col=1:2,pch='*',cex=2.5)
legend('topleft',title='TYPE',c('Normal','ALF'),pch = 21,col=1,pt.bg  = c('red','blue'))
   
   #from all
pca.dat<- prcomp(t(dat),center = T,scale. = T,cor=F)$x[,1:2]
kc <- kmeans(pca.dat, centers=2,iter.max = 20)
        #used to test the difference between the classification and given tag 
ann <- data.frame(kc$cluster,tag=c(rep('1',10),rep('2',17)))
ann$kc.cluster==ann$tag
   #plot
plot(pca.dat, col = kc$cluster,cex=1,main='K-means Clustering PCA Scatter Plot \n k=2',lwd=1.5)
points(kc$centers,col=1:2,pch='*',cex=2.5)
legend('topleft',title='TYPE',c('ALF','Normal'),pch = 21,col=1,pt.bg  = c('red','blue'))


#train
ctr <- cbind(dat.selected[,1:7],dat.selected[,11:22])
ctr <- as.data.frame(t(cbind(dat.selected[,1:7],dat.selected[,11:22])))
cte <- cbind(dat.selected[,8:10],dat.selected[,23:27])
cte <- as.data.frame(t(cbind(dat.selected[,8:10],dat.selected[,23:27])))

train.set <- as.data.frame(cbind(c(rep('Normal',7),rep('ALF',12)),ctr))
test.set <- as.data.frame(cbind(c(rep('Normal',3),rep('ALF',5)),cte))

actual.sample <- test.set[,1]
test.set <- test.set[,-1]


dat.lda <- lda(train.set[,1]~.,train.set[,2:3199])
dat.pred <- predict(dat.lda,test.set)

table(dat.pred$class,actual.sample)

plot(dat.pred$x,bg=as.numeric(factor(actual.sample)),pch=21,col=1,ylab="Discriminant 
function",axes=T,xlab="Score",main="Discriminant function for GSE38941")
legend('topright',title='TYPE',c('Normal','ALF'),pch = 21,col=1,pt.bg  = c('red','black'))                                   


# Differential expression 
normal.mean <- apply(log2(dat.selected[,1:10]),1,mean,na.rm=T)
ALF.mean <- apply(log2(dat.selected[,11:27]),1,mean,na.rm=T)
fold.value <- ALF.mean-normal.mean
fold.m <- data.frame(fold.value=fold.value,gene=row.names(dat.selected))

min(fold)

p.trans <-  -1 * log10(dat.selected.p)
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',main='Volcano Plot\nGC and ACT group 
differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(4))],fold[(p.trans> -log10(.05)&fold>log2(4))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(4))],fold[(p.trans> -log10(.05)&fold< -log2(4))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(4))
abline(h=log2(4))

max.fold <- sort(fold.m[,2],decreasing = T)
min.fold <- sort(fold.m[,2],decreasing = F)
max.fold[1:5]
max.fold[1:5]


