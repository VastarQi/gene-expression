BiocManager::install('marray')
library(dplyr)
library(limma)
library(dplyr)
library(marray)
BiocManager::install("arrayQuality")
library(arrayQuality)
data("swirl")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

datadir <- system.file("Samsfish", package = "marray")

GSM304445 <- readGPR ( fnames = "GSM304445.gpr" ,  path = 'C:/Users/97481/OneDrive/JHU/gene expression/sourcedata' )
GSM304446 <- readGPR ( fnames = "GSM304446.gpr" ,  path = 'C:/Users/97481/OneDrive/JHU/gene expression/sourcedata' )
GSM304447 <- readGPR ( fnames = "GSM304447.gpr" ,  path = 'C:/Users/97481/OneDrive/JHU/gene expression/sourcedata' )
GSM304448 <- readGPR ( fnames = "GSM304448.gpr" ,  path = 'C:/Users/97481/OneDrive/JHU/gene expression/sourcedata' )

fbGvR45 <- data.frame(Cy3MF45=GSM304445$GfMedian, Cy5MF45=GSM304445$RfMedian,Cy3MB45=GSM304445$GbMedian,Cy5MB45=GSM304445$RbMedian)
fbGvR46 <- data.frame(Cy3MF45=GSM304446$GfMedian, Cy5MF46=GSM304446$RfMedian,Cy3MB46=GSM304446$GbMedian,Cy5MB46=GSM304446$RbMedian)
fbGvR47 <- data.frame(Cy3MF47=GSM304447$GfMedian, Cy5MF47=GSM304447$RfMedian,Cy3MB47=GSM304447$GbMedian,Cy5MB47=GSM304447$RbMedian)
fbGvR48 <- data.frame(Cy3MF48=GSM304448$GfMedian, Cy5MF48=GSM304448$RfMedian,Cy3MB48=GSM304448$GbMedian,Cy5MB48=GSM304448$RbMedian)

maNorm(fbGvR45,norm = 'median')

gprfile <- read.GenePix (path='C:/study/JHU',skip=33)
gprfile <- read.GenePix (path='C:/Users/97481/OneDrive/JHU/gene expression/sourcedata/GSM/',skip=33)

gprnormMedian <- maNorm(gprfile[,1:4],norm = 'median',span=0.45)
gprnormLoess <- maNorm(gprfile[,1:4],norm = 'loess',span=0.45)
gprnormPTGloess <- maNorm(gprfile[,1:4],norm = 'printTipLoess',span=0.45)
gprnormNon <- maNorm(gprfile[,1:4],norm = 'none',span=0.45)

par(mfrow=c(4,1))
maPlot(gprnormNon,main='None Normalization')

maPlot(gprnormMedian,main='global median location normalization')
maPlot(gprnormLoess,main='loess Normalization')
maPlot(gprnormPTGloess,main='print-tip-group intensity dependent location normalization')

gprnorm04Median <- maNorm(gprfile[,4],norm = 'median',span=0.45)
gprnorm04Loess <- maNorm(gprfile[,4],norm = 'loess',span=0.45)
gprnorm04PTGloess <- maNorm(gprfile[,4],norm = 'printTipLoess',span=0.45)
gprnorm04Non <- maNorm(gprfile[,4],norm = 'none',span=0.45)
gprnorm04Median<- rm.na(gprnorm04Median)

par(mfrow=c(1,1))
plot(density(rm.na(maM(gprnorm04PTGloess[,1]))),lwd=2, col='red',
main="Density plots of log-ratios of 4 normalization in array#04",
xlab = "Standard Deviation (for data on the log2 scale)",
ylab = "Density")
lines(density(rm.na(maM(gprnorm04Median[,1]))),lwd=2,col='blue')
lines(density(rm.na(maM(gprnorm04Loess[,1]))),lwd=2,col='green')
lines(density(rm.na(maM(gprfile[,4]))),col='orange')
legend("topright",
legend = c("print-tip-group intensity","global median","loess","Pre-Normalization"),
col=c('red','blue','green','orange'),lty=1,lwd=2)


redB <- gprfile@maRb[,1:4]
redF <- gprfile@maRf[,1:4]
subRed <- redF-redB
subRed[subRed<0]=NA
logred <- log2(subRed)
redmedian<- apply(logred,2,median,na.rm=T)
redNor <- sweep(logred,2,redmedian)
2**median(redNor[,1],na.rm=T)
2**median(redNor[,2],na.rm=T)
2**median(redNor[,3],na.rm=T)
2**median(redNor[,4],na.rm=T)


colnames(redNor)<- c("array1","array2","array3","array4")
cor(redNor,method='spearman',use='complete.obs')
lossM <- gprnormLoess@maM
colnames(lossM)<- c("array1","array2","array3","array4")
cor(lossM,method='spearman',use='complete.obs')
colnames(subRed)<- c("array1","array2","array3","array4")

pairs(~array1+array2+array3+array4,panel = panel.smooth,data=redNor,main="scatter plot matrix for global median normalization")
pairs(~array1+array2+array3+array4,panel = panel.smooth,data=subRed,main="scatter plot matrix for loess normalized M values")


redB <- gprfile@maRb[,1:4]
redF <- gprfile@maRf[,1:4]
subRed <- redF-redB
subRed[subRed<0]=NA
sorted.sub.red <- apply(subRed,2,sort)
sorted.sub.red <- matrix(c(sorted.sub.red$array1,sorted.sub.red$array2,sorted.sub.red$array3,sorted.sub.red$array4),ncol = 4)
#sorted.sub.red <- cbind(sorted.sub.red$array1,sorted.sub.red$array2,sorted.sub.red$array3,sorted.sub.red$array4)
#sorted.sub.red <- matrix(sorted.sub.red)
colnames(sorted.sub.red)<- c("array1","array2","array3","array4")
sorted.sub.red.mean <- rowMeans(sorted.sub.red)

red.row.mean <- data.matrix(sorted.sub.red.mean)
red.row.mean <- data.matrix()
red.row.mean <- matrix('array1'=sorted.sub.red.mean,'array2'=sorted.sub.red.mean,'array3'=sorted.sub.red.mean,'array4'=sorted.sub.red.mean)

red.row.mean <- matrix(data = c(sorted.sub.red.mean,sorted.sub.red.mean,sorted.sub.red.mean,sorted.sub.red.mean), ncol = 4)

sub.red.rank <- apply(subRed,2,rank,ties.method='first')

array1.rank <- matrix(c(red.row.mean[,1],sub.red.rank[,1]),ncol=2 )
array2.rank <- matrix(c(red.row.mean[,2],sub.red.rank[,2]),ncol=2 )
array3.rank <- matrix(c(red.row.mean[,3],sub.red.rank[,3]),ncol=2 )
array4.rank <- matrix(c(red.row.mean[,4],sub.red.rank[,4]),ncol=2 )

array1.rank <- array1.rank[order(array1.rank[,2]),]
array2.rank <- array2.rank[order(array2.rank[,2]),]
array3.rank <- array3.rank[order(array3.rank[,2]),]
array4.rank <- array4.rank[order(array4.rank[,2]),]

ranked.mean.matrix <- cbind(array1.rank[,1],array2.rank[,1],array3.rank[,1],array4.rank[,1])
colnames(ranked.mean.matrix)<- c("array1","array2","array3","array4")
ranked.mean.matrix.log <- log2(ranked.mean.matrix)
cor(ranked.mean.matrix.log,method='spearman',use='complete.obs')
pairs(~array1+array2+array3+array4,panel = panel.smooth,data=ranked.mean.matrix.log,main="scatter plot matrix for quantile normalized ")

