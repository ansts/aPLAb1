require(stringdist)
require(parallel)
require(matrixStats)
require(igraph)
require(reshape2)
require(stringi)
require(gplots)
require(chemometrics)
require(robustbase)
require(pdist)
require(qualV)
require(msa)
require(uwot)
require(vioplot)
require(corrplot)
require(Biostrings)
require(seqinr)
require(rgl)
require(pbapply)
require(limma)

cpl3=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))
cpl1=colorRampPalette(c("#0000FF80","#00FF0080","#FFFF0080","#FF000080"))
aa=AA_ALPHABET[1:20]


# Reads from 2 sequencing experiments on phage libraries selected on control and 
# APLS serum IgM

aPL_F=read.table(file="uniqueP_QF_aLP.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
aPL_R=read.table(file="uniqueP_QR_aLP.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
aPL1_F=read.table(file="uniqueP_QF_aLP1.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
aPL1_R=read.table(file="uniqueP_QR_aLP1.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
cntr_F=read.table(file="uniqueP_QF_c1.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
cntr_R=read.table(file="uniqueP_QR_c1.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
cntr1_F=read.table(file="uniqueP_QF_c2.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)
cntr1_R=read.table(file="uniqueP_QR_c2.txt", header = FALSE,  quote=NULL, stringsAsFactors = F)

# previously published global IgM repertoire selected mimotopes
mims=read.table(file="mims.txt", header = F, quote = NULL, stringsAsFactors = F)[,1]

aPL=rbind(aPL_F,aPL_R,aPL1_F,aPL1_R)
xa=aggregate(aPL[,2], by=list(as.character(aPL[,1])), FUN=sum)
xa1=xa[xa[,2]==1,1]
xa2=xa[xa[,2]==2,1]
xa3up=xa[xa[,2]>2&xa[,2]<101,1]

cntr=rbind(cntr_F,cntr_R,cntr1_F,cntr1_R)
xc=aggregate(cntr[,2], by=list(as.character(cntr[,1])), FUN=sum)
xc1=xc[xc[,2]==1,1]
xc2=xc[xc[,2]==2,1]
xc3up=xc[xc[,2]>2&xc[,2]<101,1]

# cross-multiplicates expanding the library by including reads that are found 
# repeatedly in different experiments although in 1 or two copies per experment

xac12=intersect(xa1,xc2)
xac21=intersect(xa2,xc1)
xac123=intersect(union(xa1,xa2),xc3up)
xac321=intersect(union(xc1,xc2),xa3up)
xam123=intersect(union(xa1,xa2),mims)
xcm123=intersect(union(xc1,xc2),mims)

# the ultimate collections
aPL=unique(c(xa3up,xac12,xac21,xac123,xam123))
cntr=unique(c(xc3up,xac12,xac21,xac321,xcm123))

rm(yi,xi,aPL_F,aPL_R,aPL1_F,aPL1_R,cntr_F,cntr_R,cntr1_F,cntr1_R,xac12,xac21, xac123, xac321, xam123, xcm123)

write.csv(aPL, file="aPL.csv")
write.csv(cntr, file="cntr.csv")

# Removing Target Unrelated Peptides determined by the SAROTUP software

aPLTUPS=read.table(file="aPLTuP.txt", header = T, sep="\t",quote = NULL, stringsAsFactors = F)
excpt=unique(aPLTUPS[,3])
aPLbad=aPLTUPS[!(aPLTUPS[,3] %in% excpt[c(4,5,8)]),1]
aPL=setdiff(aPL,aPLbad)

cntrTUPS=read.table(file="cntrTuP.txt", header = T, sep="\t",quote = NULL, stringsAsFactors = F)
excpt=unique(cntrTUPS[,3])
cntrbad=cntrTUPS[!(cntrTUPS[,3] %in% excpt[c(3,5)]),1]
cntr=setdiff(cntr,cntrbad)

mimsTUPS=read.table(file="mimsTUPs.txt", header = T, sep="\t",quote = NULL, stringsAsFactors = F)
excpt=unique(mimsTUPS[,3])
mimsbad=mimsTUPS[!(mimsTUPS[,3] %in% excpt[c(3,5,6)]),1]
mims=setdiff(mims,mimsbad)
mimcl=read.csv(file="ODX79030.csv", stringsAsFactors = F)
x=mimcl[,3]
names(x)=mimcl[,2]
mimcl=x+1

Na=length(aPL)
Nc=length(cntr)
Nm=length(mims)
Nac=length(intersect(aPL,cntr))
Nam=length(intersect(aPL,mims))
Ncm=length(intersect(cntr,mims))
Nacm=length(intersect((intersect(aPL,cntr)),mims))

# Pooled library and a flag for the initial libraries each sequence appeared in
mix=unique(c(cntr, aPL, mims))
fl=1*(mix %in% cntr)+2*(mix %in% aPL)+4*(mix %in% mims)
names(fl)=mix

# Adjacency list linking sequences with longest common subsequence length > 4 
proct=proc.time()
AdjLcs=adjL(mix, t=5)
print(proc.time()-proct)
save(AdjLcs, file="AdjLcs_5")
load("Gmix")     # from the Graph

AdjLneat=AdjLcs[fl!=4] #Adjacency list of only the APLS experiment
rm(AdjLcs)
cl=makeCluster(4)
clusterExport(cl, c("fl"))
AdjLneat=pblapply(AdjLneat, function(l){
  l[fl[l]!=4]
}, cl=cl)
stopCluster(cl)
AdjLneat=AdjLneat[lengths(AdjLneat)>0]
save(AdjLneat, file="AdjLneat")
rm(AdjLneat)

t0=rep(0,7)
names(t0)=1:7

# Building a table of neighbor distributions for each vertex across diagnosis flag
proct=proc.time()
cl = makeCluster(getOption("cl.cores", 4))
clusterExport(cl,varlist=c("fl","mix","t0"), envir = environment())
Ngblcs=t(parSapply(cl, AdjLcs, function(v){
  x=table(fl[mix %in% v])
  t=t0
  t[names(x)]=x
  return(t)
}))
stopCluster(cl)
print(proc.time()-proct)

# Regrouping
n=names(AdjLneat)
Ngbln=Ngblcs[n,]
Ngbln[,1]=round((Ngbln[,1]+Ngbln[,5]+Ngbln[,7]+(Ngbln[,3]/2)))
Ngbln[,2]=round((Ngbln[,2]+Ngbln[,6]+Ngbln[,7]+(Ngbln[,3]/2)))
Ngbln=Ngbln[,-c(3,4,5,6,7)]
Ngbln=Ngbln[rowSums(Ngbln)>0,]
fln=fl[rownames(Ngbln)]
fln[fln==5]=1
fln[fln==6]=2
fln[fln==7]=3
tfln=table(fln)
pnn=(tfln[1:2]+tfln[3]/2)/sum(tfln) # the neighboring sequences which are found 
                                    # in both libraries should count to both equally
                                    # but the degree of the vertex is preserved 
                                    # if the counts are split in half

cl=makeCluster(4)
clusterExport(cl, "pnn")
pNNNn=pbapply(Ngbln,1,function(l){
  s=sum(l)
  if (s>30) prop.test(l[1], s, p=pnn[1])$p.value else binom.test(l[1], s, p=pnn[1])$p.value
}, cl=cl)
stopCluster(cl)
adj.pNNNn=p.adjust(pNNNn)
# vulcano plot of contr vs aPL
pdf(file="Vulcano_aPLvsContr_newn.pdf", height=10,width = 10)
plot(log(((Ngbln[,2]+0.5)*pnn[1])/((Ngbln[,1]+0.5)*pnn[2])), -log(pNNNn), pch=16, col=1+(adj.pNNNn<0.05)+(log(Ngbln[,2]/Ngbln[,1])<1&adj.pNNNn<0.05)*1, cex=-log(pNNNn)/45+0.01, xlab="log Fold Predominance", ylab="-log(p)", xlim=c(-3,3))
dev.off()

aPLScalls=Ngbln[(((Ngbln[,2]+0.5)*pnn[1])/((Ngbln[,1]+0.5)*pnn[2])>1)&adj.pNNNn<0.05,]
aPLScalls=cbind(aPLScalls, (aPLScalls[,2]*pnn[1])/(aPLScalls[,1]*pnn[2]), adj.pNNNn[rownames(aPLScalls)])
aPLScalls=aPLScalls[order(aPLScalls[,3],decreasing = T),]
CntrCalls=Ngbln[(((Ngbln[,2]+0.5)*pnn[1])/((Ngbln[,1]+0.5)*pnn[2])<1)&adj.pNNNn<0.05,]
CntrCalls=cbind(CntrCalls, (CntrCalls[,2]*pnn[1])/(CntrCalls[,1]*pnn[2]), adj.pNNNn[rownames(CntrCalls)])
CntrCalls=CntrCalls[order(CntrCalls[,3],decreasing = T),]
save(aPLScalls,CntrCalls, file="vulcalls")

yr=range(-log(pNNNn))
gr=fln[adj.pNNNn<0.05]
y=log((Ngbln[adj.pNNNn<0.05,2]*pnn[1])/(Ngbln[adj.pNNNn<0.05,1]*pnn[2]))
grr=c("Contr.","aPL","Contr.+aPL")
grrs=grr[gr]
vioplot(y[grrs=="Contr."],y[grrs=="aPL"],y[grrs=="Contr.+aPL"], names=c("Contr.","aPL","Contr.+aPL"), horizontal = T, ylim=c(-3,3))

