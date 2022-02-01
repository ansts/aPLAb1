i=commandArgs(trailingOnly = TRUE)
load(file="varsNgbTGlo") # varsNgbTI
source("ngbtbl.R")
i=as.numeric(i)
i1=(i-1)*16+1
i2=i*16

j=seq_along(AdjL_Glo) # AdjL_I
ji=cut(j, 160, labels = F)
j=which(ji %in% i1:i2)
NgbTblGlo=ngbtbl(AdjL_Glo[j],  fl)  #NgbTblI, AdjL_I,flI
save(NgbTblGlo, file=paste("NgbTblGlo","part",i,sep="_",collapse="")) # NgbTblI, NgbTblI