rat_kd <- read.table('C:/Users/97481/OneDrive/JHU/gene expression/sourcedata/rat_KD.txt',header=T,row.names=1)
rat_kd <- log2(rat_kd)
names(rat_kd)
control <- rat_kd[,(1:6)]
keto <- rat_kd[,(7:11)]

t.test.all.genes <- function(x,s1,s2) {
  x1 <- as.vector(s1)
  x2 <- as.vector(s2)
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

rat_kd[control]

as.vector(control)
t.test.all.genes(rat_kd,control,keto)
pv <- apply(rat_kd,1,t.test.all.genes,s1=control,s2=keto)


t.test.all.genes <- function(x,s1,s2) {
  x1 <- as.numeric(unlist(s1))
  x2 <- as.numeric(unlist(s2))
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}

control <- as.character(names(rat_kd[,(1:6)]))
keto <- as.character(names(rat_kd[,(7:11)]))