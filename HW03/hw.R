BiocManager::install('multtest')
library(multtest)
library(Biobase)
library(annotate)
library(limma)
data(golub)

golub.datfr <- as.data.frame(golub)
growname <- paste('g',1:nrow(golub.datfr),sep = '')
rownames(golub.datfr) <- growname
colnames(golub.datfr) <- paste(colnames(golub.datfr),golub.cl,sep = ' ')

wilcox.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- wilcox.test(x1,x2,exact=F,alternative="two.sided",correct=T)
  out <- as.numeric(t.out$statistic)
  return(out)
}


original.wmw.run <- apply(golub.datfr,1,wilcox.test.all.genes,s1=golub.cl==0,s2=golub.cl==1)

wmw.max.list <- c()
for (i in 1:500) {
  golub.datfr <- golub.datfr[,sample(ncol(golub.datfr))]
  wmw.process <- apply(golub.datfr, 1, wilcox.test.all.genes, s1=c(1:27), s2=c(28:38))
  wmw.max.list<-c(wmw.max.list, max(wmw.process))
}

original.wmw.run.95.max <- original.wmw.run[original.wmw.run >quantile(wmw.max.list, probs=0.95)]
summary(original.wmw.run.95.max)
attributes(original.wmw.run.95.max)
summary(original.wmw.run.95.max)

design <-cbind(Grp1=1,Grp2vs1=c(rep(0,length(golub.cl[golub.cl==0])),rep(1,length(golub.cl[golub.cl==1]))))
fit <- lmFit(golub.datfr,design)
eb.fit <- eBayes(fit)
attributes(eb.fit)
eb.p.value <- eb.fit$p.value[,2]

s.eb.p.value <- sort(eb.p.value)[1:length(original.wmw.run.95.max)]
intersect(names(s.eb.p.value),names(original.wmw.run.95.max))

t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2,alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

t.result <- apply(golub.datfr, 1, t.test.all.genes, s1=golub.cl==0, s2=golub.cl==1)
t.result.95.max <- t.result[t.result < 0.01]

plot(t.result.95.max, eb.p.value[names(t.result.95.max)], xlab='P-Values of Students T-Test ',ylab='P-Values of Empirical Bayes Method', 
     main='P-Value Comparison Plot\n StudentsT-Test vs. Empirical Bayes \nGolub Data(P-Value<= 0.01)', col='black', bg='red', pch=21, cex=0.5)





