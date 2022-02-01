i=commandArgs(trailingOnly = TRUE)
load(file="varsAdja") #varsAdj, varsAdjI
source("adjL.R")
i=as.numeric(i)
i1=(i-1)*16+1
i2=i*16

#s=which(flI!=4)  #fl

j=seq_along(fla) #s
ji=cut(j, 160, labels = F)
j=which(ji %in% i1:i2)
AdjLcs=adjL(mixa, t=5, subs=j) #mix, s[j]
save(AdjLcs, file=paste("AdjL_a","part",i,sep="_",collapse="")) # AdjLcs
