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

index_to_mean <- function(input.rank, input.mean){
  return(input.mean[input.rank])
}

q.norm.result <- apply(sub.red.rank, 2, index_to_mean, input.mean=sorted.sub.red.mean)
q.norm.result

par(mfrow=c(4,1))
hist(q.norm.result[,1])
hist(q.norm.result[,2])
hist(q.norm.result[,3])
hist(q.norm.result[,4])

q.norm.result.log <- log2(q.norm.result)
cor(q.norm.result.log,method='spearman',use='complete.obs')
pairs(~array1+array2+array3+array4,panel = panel.smooth,data=q.norm.result.log,main="scatter plot matrix for quantile normalized ")
