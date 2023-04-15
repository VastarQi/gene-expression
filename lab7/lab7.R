install.packages('multtest')
BiocManager::install('multtest')
BiocManager::install('ggfortify')
library(ggfortify)
library(MASS)
library
#2
Sotiriou  <- read.table('C:/study/JHU/Sotiriou.txt',header=T,row.names=1)
Sotiriou.ann<- read.table('C:/study/JHU/Sotiriou_annotations.txt',header=T,row.names=1)
Sotiriou.pca <- prcomp(t(Sotiriou),cor=F)
Sotiriou.pca.loading <- Sotiriou.pca$x[,1:2]

plot(range(Sotiriou.pca.loading[,1]),range(Sotiriou.pca.loading[,2]),type="n",xlab='p1',ylab='p2',
     main='PCA plot of Sotiriou breast cancer data set\np2 vs. p1')
points(Sotiriou.pca.loading[,1][Sotiriou.ann$site=='KIU'], Sotiriou.pca.loading[,2][Sotiriou.ann$site=='KIU'],
       col=1,bg='red',pch=21,cex=1)
points(Sotiriou.pca.loading[,1][Sotiriou.ann$site=="OXF"], Sotiriou.pca.loading[,2][Sotiriou.ann$site=="OXF"],
       col=1,bg='blue',pch=21,cex=1)
legend('topright',title='site',c('KIU','OXF'),pch = 21,col=1,pt.bg  = c('red','blue'))




#3
Sotiriou.pca.var <- round(Sotiriou.pca$sdev^2 / sum(Sotiriou.pca$sdev^2)*100,2)
plot(c(1:length(Sotiriou.pca.var)),Sotiriou.pca.var,type="b",
     xlab="# components",ylab="% variance",pch=21,col=1,bg=3,cex=1)
title("Scree plot showing % variability explained by each eigenvalue\nSotiriou breast cancer dataset")

#4
#non-metric MDS
Sotiriou.dist <- dist(t(Sotiriou))
Sotiriou.mds <- isoMDS(Sotiriou.dist)
plot(Sotiriou.mds$points, type = "n",xlab='p1',ylab='p2')
points(Sotiriou.mds$points[,1][Sotiriou.ann$site=='KIU'], Sotiriou.mds$points[,2][Sotiriou.ann$site=='KIU'],
       col='red',pch=16,cex=1)
points(Sotiriou.mds$points[,1][Sotiriou.ann$site=="OXF"], Sotiriou.mds$points[,2][Sotiriou.ann$site=="OXF"],
       col='blue',pch=16,cex=1)
title(main='non-metric MDS plot of Sotiriou breast cancer data set\nstress=20%')
legend('bottomright',title='site',c('KIU','OXF'),col=1,pt.bg  = c('red','blue'),pch=21,cex=1,horiz=F)

#classic
Sotiriou.loc <- cmdscale(Sotiriou.dist)
plot(Sotiriou.loc, type = "n",xlab='p1',ylab='p2')
points(Sotiriou.loc[,1][Sotiriou.ann$site=='KIU'], Sotiriou.loc[,2][Sotiriou.ann$site=='KIU'],
       col='red',pch=16,cex=1)
points(Sotiriou.loc[,1][Sotiriou.ann$site=="OXF"], Sotiriou.loc[,2][Sotiriou.ann$site=="OXF"],
       col='blue',pch=16,cex=1)
title(main='classic MDS plot of Sotiriou breast cancer data set')
legend('bottomright',title='site',c('KIU','OXF'),col=1,pt.bg  = c('red','blue'),pch=21,cex=1,horiz=F)

#5
temp <- t(Sotiriou)
temp <- scale(temp,center=T,scale=T) 

k.speClust2 <- function (X, qnt=NULL) {
  dist2full <- function(dis) {
    n <- attr(dis, "Size")
    full <- matrix(0, n, n)
    full[lower.tri(full)] <- dis
    full + t(full)
  }
  dat.dis <- dist(t(X),"euc")^2
  if(!is.null(qnt)) {eps <- as.numeric(quantile(dat.dis,qnt))}
  if(is.null(qnt)) {eps <- min(dat.dis[dat.dis!=0])}
  kernal <- exp(-1 * dat.dis/(eps))
  K1 <- dist2full(kernal)
  diag(K1) <- 0
  D = matrix(0,ncol=ncol(K1),nrow=ncol(K1))
  tmpe <- apply(K1,1,sum)
  tmpe[tmpe>0] <- 1/sqrt(tmpe[tmpe>0])
  tmpe[tmpe<0] <- 0
  diag(D) <- tmpe
  L <- D%*% K1 %*% D
  X <- svd(L)$u
  Y <- X / sqrt(apply(X^2,1,sum))
}
phi <- k.speClust2(t(temp),qnt=NULL)
plot(range(phi[,1]),range(phi[,2]),xlab="phi1",ylab="phi2",
     main="Weighted Graph Laplacian plot\nSotiriou breast cancer dataset")
points(phi[,1][Sotiriou.ann$site=='KIU'],phi[,2][Sotiriou.ann$site=='KIU'],col="red",pch=16,cex=1)
points(phi[,1][Sotiriou.ann$site=="OXF"],phi[,2][Sotiriou.ann$site=="OXF"],col="blue",pch=16,cex=1)
legend('left',title='site',c("KIU", "OXF"),col=c("red", "blue"),pch=16,cex=.7,horiz=F)


#test
autoplot(Sotiriou.pca)
plot(range(Sotiriou.pca$x[,1]),range(dat.loadings[,2]),type="n",xlab='p1',ylab='p2',main='PCA plot of Golub Data\np2 vs. p1')
points(dat.loadings[,1][ann.dat2==1], dat.loadings[,2][ann.dat2==1],col=1,
       bg='red',pch=21,cex=1.5)
points(dat.loadings[,1][ann.dat2==0], dat.loadings[,2][ann.dat2==0],col=1,bg='blue',pch=21,cex=1.5)

Sotiriou.pca.loading[,2][Sotiriou.ann$site=='KIU'] <- "red"
Sotiriou.pca.loading[,2][Sotiriou.ann$site=="OXF"] <- "blue"

plot(Sotiriou.pca,col=Sotiriou.ann,xlab="p1",ylab="p2",pch=16,cex=1.5,main="Kernel PCA of Golub data\nsigma=0.002")
legend(1.5,1.5,c("AML", "ALL"),col=c("red", "blue"),pch=15,cex=.7,horiz=F)
