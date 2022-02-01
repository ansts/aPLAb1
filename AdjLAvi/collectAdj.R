require(parallel)
require(pbapply)

load("AdjLcs_part_1")
AdjL_Glo=AdjLcs
for (i in 2:10){
  load(paste("AdjLcs_part_",i,sep=""))
  AdjL_Glo=c(AdjL_Glo,AdjLcs)
}
save(AdjL_Glo,file="AdjL_Glo")
load("fl")
fls=fl[fl!=4]
load("AdjLcsI_part_1")
AdjL_I=AdjLcs
for (i in 2:10){
  load(paste("AdjLcsI_part_",i,sep=""))
  AdjL_I=c(AdjL_I,AdjLcs)
}
save(AdjL_I,file="AdjL_I")

load("flI")
load("AdjL_I")
save(flI,AdjL_I, file="varsNgbTI")


#  this part was sent to HPC cluster --------------------------------------


# t0=rep(0,7)
# names(t0)=1:7
# proct=proc.time()
# cl=makeCluster(4)
# clusterExport(cl, c("flI", "t0"), envir = environment())
# NgbI=t(pbsapply(AdjL_I, function(l){
#   lf=flI[l]
#   t1=table(lf)
#   t0[names(t1)]=t1
#   return(t0)
# }, cl=cl))
# stopCluster(cl)
# print(proc.time()-proct)
# save(NgbI, file="NgbI")


# -------------------------------------------------------------------------


rm(AdjL_I, flI)
load("AdjL_Glo")
load("fl")
save(fl,AdjL_Glo, file="varsNgbTGlo")


#  this part was sent to HPC cluster --------------------------------------



# proct=proc.time()
# cl=makeCluster(4)
# clusterExport(cl, c("fln", "t0"), envir = environment())
# NgbGlo=t(pbsapply(AdjL_Glo, function(l){
#   lf=fln[l]
#   t1=table(lf)
#   t0[names(t1)]=t1
#   return(t0)
# }, cl=cl))
# stopCluster(cl)
# print(proc.time()-proct)
# save(NgbGlo, file="NgbGlo")
# rm(fln,AdjL_Glo)

load("NgbTblGlo_part_1")
NgbG=NgbTblGlo
for (i in 2:10){
  load(paste("NgbTblGlo_part_",i,sep=""))
  NgbG=rbind(NgbG,NgbTblGlo)
}
save(NgbG,file="NgbG")

load("NgbTblI_part_1")
NgbI=NgbTblI
for (i in 2:10){
  load(paste("NgbTblI_part_",i,sep=""))
  NgbI=rbind(NgbI,NgbTblI)
}
save(NgbI,file="NgbI")

x=NgbI
x[,1]=x[,1]+x[,5]
x[,2]=x[,2]+x[,6]
x[,3]=x[,3]+x[,7]
y=x
all(x[,1:3]==y[,1:3])
Ngba=x[,1:3]
save(Ngba,file="Ngba")

load("AdjL_a_part_1")
AdjL_a=AdjLcs
for (i in 2:10){
  load(paste("AdjL_a_part_",i,sep=""))
  AdjL_a=c(AdjL_a,AdjLcs)
}
save(AdjL_a,file="AdjL_a")

     