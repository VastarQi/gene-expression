library(MASS)

dat <- read.table('C:/Users/QiHY//OneDrive/JHU/gene expression/sourcedata/lung_cancer_data.txt', header = T , row.names = 1)
classname <- c(rep('Adeno',10),rep('SCLC',9),rep('Normal',5))
dat.t <- data.frame(classname,t(dat))

train.set <- rbind(dat.t[1:6,],dat.t[11:16,],dat.t[20:22,])
test.set <- rbind(dat.t[7:10,],dat.t[17:19,],dat.t[23:24,])
actual.sample <- test.set[,1]
test.set <- test.set[,-1]


dat.lda <- lda(train.set[,1]~.,train.set[,2:3])
dat.pred <- predict(dat.lda,test.set[,1:2])
table(dat.pred$class,actual.sample)

plot(dat.pred$x,bg=as.numeric(factor(actual.sample)),pch=21,col=1,ylab="Discriminant 
function",axes=T,xlab="Score",main="Discriminant function for Colon dataset")
legend()                                   

