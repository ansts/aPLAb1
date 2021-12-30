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
require(emmeans)
require(Rcapture)
require(eulerr)

cpl3=colorRampPalette(c("#AFAFAF0A","#FF0FFF0A","#00FFF00A"))
cpl1=colorRampPalette(c("#0000FF80","#00FF0080","#FFFF0080","#FF000080"))
cpl0=colorRampPalette(c("#0000FF80","#0000FF80","#0000FF80","#00008080","#B0FF0080","#FFFF0080","#FF000080"))
aa=AA_ALPHABET[1:20]


# Raw data from sequencing ------------------------------------------------

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

# previously published (old) and a second unpublished experiment 
# derived global IgM repertoire selected mimotopes, cleaning TUPs
IgM1=load("oldmims")
IgM2=load("IVIgMseqbeforeTUPclean")
allIgMTUPs=read.table("exp1.txt",sep = '\t', header = T)
oldmimsTUPs=read.table("oldmimsTUPs.txt",sep = '\t', header = T)
typbad=unique(c(allIgMTUPs$Brief.Description,oldmimsTUPs$Brief.Description))
aIbad=allIgMTUPs[allIgMTUPs$Brief.Description %in% typbad[c(1,3,4,5,7)], 1]
oldbad=oldmimsTUPs[oldmimsTUPs$Brief.Description %in% typbad[c(1,3,4,5,7)], 1]
aIgM=IVIgMseq[!(IVIgMseq$Seq %in% aIbad),]
oldm=oldmims[!(oldmims$Seq %in% oldbad),]
allseq=rbind(aIgM,oldm)
allseq=aggregate(allseq$N, by=list(allseq$Seq), "sum")
colnames(allseq)=c("Seq","N")
allseq=allseq[order(allseq$N, decreasing=T),]
mims=allseq$Seq[allseq$N>2&allseq$N<101]
rm(allseq, oldm, aIgM, aIbad,oldbad, typbad, IgM1, IgM2)

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
# repeatedly in different experiments although in 1 or two copies per experiment

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

Na=length(aPL)
Nc=length(cntr)
Nm=length(mims)
Nac=length(intersect(aPL,cntr))
Nam=length(intersect(aPL,mims))
Ncm=length(intersect(cntr,mims))
Nacm=length(intersect((intersect(aPL,cntr)),mims))

# Pooled library and a flag for the initial libraries each sequence appeared in
mix=unique(c(cntr, aPL, mims))
mixa=unique(c(cntr, aPL))
fl=1*(mix %in% cntr)+2*(mix %in% aPL)+4*(mix %in% mims)
fla=1*(mixa %in% cntr)+2*(mixa %in% aPL)
names(fl)=mix
names(fla)=mixa


# Estimates of original diversity -----------------------------------------

# From:
# W. L. Matochko, R. Derda. "Error analysis of deep sequencing of phage libraries: peptides censored in 
# sequencing". Comput Math Methods Med, 2013, 491612. 2013.
# and
# Matochko WL, Cory Li S, Tang SK, Derda R. Prospective identification of parasitic sequences in phage 
# display screens. Nucleic Acids Res. 2014;42(3):1784-98. Epub 2013/11/13. doi: 10.1093/nar/gkt1104. 
# PhD7-Naive-30FuR.txt
#
# http://www.chem.ualberta.ca/∼derda/mathbiology/PhD7-Amp-30F.txt
#

matochkoAF30=read.table(file = "PhD7-Amp-30F.txt", stringsAsFactors = FALSE)
matochkoNF30=read.table(file = "PhD7-Naive-30FuR.txt", stringsAsFactors = FALSE)
matF30=rbind(matochkoAF30[,2:3],matochkoNF30[,2:3])
matF30=aggregate(matF30[,2], by=list(matF30[,1]), sum)
matlib=matF30[,1]
matlib1=matF30[matF30[,2]==1,1]
matF30=matF30[matF30[,2]>2&matF30[,2]<101,]
matF30=matF30[order(matF30[,2], decreasing=T),]
match7astrxi=which(!is.na(stri_locate_first_fixed(matF30[,1], "*")[,1]))
matF30=matF30[-match7astrxi,]
tups=read.table("matF30TUPs.txt", header = T, sep="\t")
typbad=unique(tups$Brief.Description)
pepbad=tups[tups$Brief.Description %in% typbad[c(1,3,4,5,7,9,11)], 1]
matpep=matF30[!matF30[,1] %in% pepbad,1]
rm(matochkoAF30, matF30,matochkoNF30,matochkoAF30,matlib,matlib1)

MvsR=intersect(mix, matpep)
GandA=sum(fl %in% 6:7)
GandC=sum(fl %in% c(5,7))
Gall=sum(fl %in% 4:7)
Aall=sum(fl %in% c(2,3,6,7))
Call=sum(fl %in% c(1,3,5,7))

LM=length(matpep)
LS=length(mix)

cRcRM=rbind(c(1,1,length(MvsR)),c(1,0,length(setdiff(matpep,mix))),c(0,1,length(setdiff(mix,matpep))))
RMmodel=closedpCI.t(cRcRM,dfreq = T,m="Mt") # Actual diversity of Ph.D.-7

tbl=table(fla)
prop.test(tbl[3],(tbl[1]+tb[3]), (tbl[2]+tbl[3])/RMmodel$results[1])
prop.test(GandC,Call, Gall/RMmodel$results[1])

tbl[3]/((tbl[1]+tb[3])/RMmodel$results[1]*(tbl[2]+tb[3]))

GandC/(Gall/RMmodel$results[1]*Call)
GandA/(Gall/RMmodel$results[1]*Aall)

aap=read.table("aaPh_ok.txt")
colnames(aap)=aap[1,]
aap=aap[2,]

cl=makeCluster(4)
clusterExport(cl,c("mix","matpep","rndpeppre"), envir = environment())
clusterEvalQ(cl, require(Rcapture))
crcmodBS0=pbsapply(1:200, function(i){  # diversity estimate based on uniform distr.
  print(i)                           # bootstrapped  
  x=rndpeppre(length(matpep))
  Is=intersect(mix, x)
  D=rbind(c(1,1,length(Is)),c(1,0,length(setdiff(x,mix))),c(0,1,length(setdiff(mix,x))))
  M=closedpCI.t(D,dfreq = T,m="Mt")
  M$CI[1]
},cl=cl)
stopCluster(cl)

summary(crcmodBS0)
hist(crcmodBS0)
quantile(crcmodBS0,c(0.05,0.95))

cl=makeCluster(4)
clusterExport(cl,c("mix","matpep","rndpeppre","aap"), envir = environment())
clusterEvalQ(cl, require(Rcapture))
crcmodBS1=pbsapply(1:200, function(i){ # diversity estimate based on PDL distr.
  print(i)                           # bootstrapped  
  x=rndpeppre(length(matpep), aa=aap)
  Is=intersect(mix, x)
  D=rbind(c(1,1,length(Is)),c(1,0,length(setdiff(x,mix))),c(0,1,length(setdiff(mix,x))))
  M=closedpCI.t(D,dfreq = T,m="Mt")
  M$CI[1]
}, cl=cl)
stopCluster(cl)

summary(crcmodBS1)
hist(crcmodBS1)
quantile(crcmodBS1,c(0.05,0.95))

cl=makeCluster(4)
clusterExport(cl,c("mix","matpep","rndpeppre"), envir = environment())
clusterEvalQ(cl, require(Rcapture))
crcmodBS2=pbsapply(1:200, function(i){ # diversity estimate based on scrambled seq
  print(i)                              # bootstrapped  
  x=apply(apply(t(sapply(strsplit(matpep, split=""), unlist)),2,sample),1,paste,sep="",collapse="")
  Is=intersect(mix, x)
  D=rbind(c(1,1,length(Is)),c(1,0,length(setdiff(x,mix))),c(0,1,length(setdiff(mix,x))))
  M=closedpCI.t(D,dfreq = T,m="Mt")
  M$CI[1]
}, cl=cl)
stopCluster(cl)

summary(crcmodBS2)
hist(crcmodBS2)
quantile(crcmodBS2,c(0.05,0.95))

# Venn

Vennall=cbind(fl %in% c(1,3,5,7), fl %in% c(2,3,6,7), fl %in% c(4:7))
Vennall=euler(Vennall, shape="ellipse", extraopt=T)
plot(Vennall, labels=c("Control", "APS", "Global"), quantities=T, adjust_labels=T)

.# The Graph ---------------------------------------------------------------

# Adjacency list linking sequences with longest common subsequence length > 4 
# only the neighborhoods of the cntr and aPL are listed for brevity

# varsAdj for Avitohol
save(fl,mix, file="varsAdj")
save(fla,mixa, file="varsAdja")
save(fl, file="fl")
save(fla, file="fla")

load("IgJt7")
IgJ=unique(IgJt7)
rm(IgJt7)
mixI=c(mixa,IgJ)
mixI=unique(mixI)
flI=1*(mixI %in% cntr)+2*(mixI %in% aPL)+4*(mixI %in% IgJ)
names(flI)=mixI
save(flI,mixI, file="varsAdjI")
save(flI, file="flI")

# from HPC Avitohol

load("NgbG") # Neighborhood table of cntr and APLS together with global repertoire IgM mimotopes
load("NgbI") # Neighborhood table of cntr and APLS together with IgJ (idiotopes)
load("Ngba") # Neighbothood table of cntr and APLS only (columns 1:3 of any of the previous + columns 5:6)
load("AdjL_a")
colnames(NgbG)=c("Cntr", "APLS", "All","Global","Global&Cntr", "Global&APLS", "Global&All")
colnames(NgbI)=c("Cntr", "APLS", "All","IgJ","IgJ&Cntr", "IgJ&APLS", "IgJ&All")
colnames(Ngba)=c("Cntr", "APLS", "All")

Ga=adjL2Glite(AdjL_a)
save(Ga, file="Ga")

# Regrouping

Ng1=Ngba
Ng1[,1]=round(Ng1[,1]+Ng1[,3]/2)
Ng1[,2]=round(Ng1[,2]+Ng1[,3]/2)
Ng1=Ng1[,-3]
colnames(Ng1)=c("Cntr","aPLS")
Ng1=Ng1[rowSums(Ng1)>0,]
tfl=table(fla) # general distribution of the cntr, aPLS and common mimos
ns=c(tfl[1]+tfl[3]/2,tfl[2]+tfl[3]/2) # regrouping corresponding to the above
pn=ns/sum(ns) # marginal p

cl=makeCluster(4)
clusterExport(cl, "pn")
pNn=pbapply(Ng1,1,function(l){
  s=sum(l)
  if (s>30) prop.test(l[1], s, p=pn[1])$p.value else binom.test(l[1], s, p=pn[1])$p.value
}, cl=cl)
stopCluster(cl)
adj.pNn=p.adjust(pNn)
x=Ng1[adj.pNn<0.05,]
hist(rowSums(Ng1), breaks=50, xlim=c(0,2000), freq = F, col=rgb(0.3,0.3,0.3,0.75), ylim=c(0,0.004), main="", xlab="Number of neighbors ")
par(new=T)
hist(rowSums(x), breaks=50, xlim=c(0,2000), col=rgb(1,0,0,0.5), freq=F, ylim=c(0,0.004), main="", xlab="", ylab="")

hist(rowSums(x)[x[,1]/pn[1]>x[,2]/pn[2]], breaks = 20, xlim=c(0,2000), col=rgb(0,0,0,0.5), freq=F, ylim=c(0,0.004), main="", xlab="", ylab="")
par(new=T)
hist(rowSums(x)[x[,1]/pn[1]<x[,2]/pn[2]], breaks = 20, xlim=c(0,2000), col=rgb(1,0,0,0.5), freq=F, ylim=c(0,0.004), main="", xlab="", ylab="")
legend("top", legend=c("Control","APS"), fill=c(rgb(0,0,0,0.5),rgb(1,0,0,0.5)), bty="n")

range(x[x[,1]/pn[1]<x[,2]/pn[2],1])
range(x[x[,1]/pn[1]>x[,2]/pn[2],1])
range(x[x[,1]/pn[1]<x[,2]/pn[2],2])
range(x[x[,1]/pn[1]>x[,2]/pn[2],2])

x=x/rowSums(x)/cbind(rep(pn[1],nrow(x)),rep(1-pn[1],nrow(x)))
xp=colSums(x>1)/sum(x>1)
# 38.3% of the significant neighborhoods contain predominantly cntr 
# while the probability of finding a cntr is pn[1]=26.1% - so?
# the significant are disproportionately more in the cntr subset?

prop.test(137, 358, p=pn[1])

# vulcano plot of contr vs aPL
pdf(file="Vulcano_aPLvsContr_u.pdf", height=10,width = 10)
plot(log(((Ng1[,2]+0.5)*pn[1])/((Ng1[,1]+0.5)*pn[2])), -log(pNn), pch=16, col=1+(adj.pNn<0.05)+(log(Ng1[,2]/Ng1[,1])<1&adj.pNn<0.05)*1, cex=-log(pNn)/45+0.01, xlab="log Fold Predominance", ylab="-log(p)", xlim=c(-3,3))
dev.off()

LOR=log10(x[,2]/x[,1])
hist(LOR, breaks=100, xlim=c(-1,1), main="Distribution of log odds of significant neighborhoods")
text(c(-0.75, 0.75), c(20, 20), labels=c("Control", "APLS"))

vulcalls=sort(LOR)
save(vulcalls, file="vulcalls")

fln=fla[rownames(Ng1)]
yr=range(-log(pNn))
gr=fln[adj.pNn<0.05]
y=log((Ng1[adj.pNn<0.05,2]*pn[1])/(Ng1[adj.pNn<0.05,1]*pn[2]))
grr=c("Contr.","aPL","Contr.+aPL")
grrs=grr[gr]
vioplot(y[grrs=="Contr."],y[grrs=="aPL"],y[grrs=="Contr.+aPL"], names=c("Contr.","aPL","Contr.+aPL"), horizontal = T, ylim=c(-3,3))

# NgbG vulcalls distr Fig

vulcallsG=(NgbG[names(LOR), 4:7]+0.1)/rowSums(Ng1[names(LOR),])
for (i in 1:4){
  plot(log10(vulcallsG[,i]), LOR, col=rgb(i %% 2, i==2, i>2,0.5), pch=16, cex=0.5, xlim=c(-3,1.5), ylim=c(-1,1), xlab="Log ratio to number of Control+APLS mimotopes in the neighborhood", main="Global specificities")
  par(new=T)
}
i=1:4
legend("bottomright", legend=colnames(NgbG)[4:7],col=rgb(i %% 2, i==2, i>2,0.5),pch=16,pt.cex = 1.5, cex=0.75)

pdf(file="vulcallsG.pdf", width = 3, height=5)
for (i in 1:4) boxplot(log10(vulcallsG[,i])~(LOR>0), names=c("Control", "APLS"), notch=T, cex=0.3, xlab="", ylab="",  xaxt='n', yaxt='n', ann=FALSE, col=rgb(i %% 2, i==2, i>2,0.5)) #main=colnames(vulcallsG)[i], ylab="Proportion of global repertoire mimotopes in the neighborhood"
dev.off()

LORplus=data.frame(LOR=LOR[LOR>0],vulcallsG[LOR>0,])
LORpL=melt(LORplus, id.vars = 1)
lmLORp=lm(LOR~0+variable*value, data=LORpL)
pairs(emtrends(lmLORp, ~variable, var="value"))
LORminus=data.frame(LOR=LOR[LOR<0],vulcallsG[LOR<0,])
LORmL=melt(LORminus, id.vars = 1)
lmLORm=lm(LOR~0+variable*value, data=LORmL)
pairs(emtrends(lmLORm, ~variable, var="value"))

# NgbI vulcalls distr Fig

vulcallsI=(NgbI[names(LOR), 4:7])/rowSums(Ng1[names(LOR),])+0.1
xr=range(log10(vulcallsI[,1:3]))
for (i in 1:3){
  plot(log10(vulcallsI[,i]), LOR, col=rgb(i %% 2, i==2, i>2,0.5), pch=16, cex=0.5, xlim=xr, ylim=c(-1,1), xlab="Log ratio to number of Control+APLS mimotopes in the neighborhood", main="Idiotopes")
  par(new=T)
}
i=1:3
legend("bottomright", legend=colnames(NgbI)[4:7],col=rgb(i %% 2, i==2, i>2,0.5),pch=16,pt.cex = 1.5, cex=0.75)

pdf(file="vulcallsI.pdf", width = 6.2, height=10)
for (i in 1:3) boxplot(log10(vulcallsI[,i])~(LOR>0), notch=T, col=rgb(i %% 2, i==2, i>2,0.5), xlab="", ylab="Log proportion of idiotopes in the neighborhood", names=c("Control", "APLS"))
dev.off()

plot(log10(NgbG[,4]),log10(NgbI[,4]))

summary(lm(log10(NgbG[,4]+0.3)~log10(NgbI[,4]+0.3)))
nG=NgbG[,5:7]
nI=NgbI[,5:7]
nGI=cbind(nG, nI)
apply(nGI,2 ,range)
nGI=nGI[rowSums(nGI)>0,]
vcn=names(vulcalls)
vcn=vcn[vcn %in% rownames(nGI)]

pdf(file="GIcors.pdGf", width=10, height=6.2)
paro=par(mfrow=c(1,2))
corrplot(cor(nGI), method="color", addCoef.col = T, cl.lim = c(0,1), col=cpl0(200))
corrplot(cor(nGI[names(vulcalls),]), method="color", addCoef.col = T, cl.lim = c(0,1), col=cpl0(200))
par(paro)
dev.off()

nGInoD=nGI/rowSums(nGI)
pdf(file="GIcorsnoD.pdf", width=10, height=6.2)
paro=par(mfrow=c(1,2))
corrplot(cor(nGInoD), method="color", addCoef.col = T)
corrplot(cor(nGInoD[names(vulcalls),]), method="color", addCoef.col = T)
par(paro)
dev.off()

IbyGlm=lm(log10(NgbI[,4]+0.5)~log10(NgbG[,4]+0.5)*log10(rowSums(Ngba+0.5)))
summary(IbyGlm)

col1=log10(rowSums(Ngba+0.5))
col1=cut(col1, 200, labels = F)
plot(log10(NgbG[,4]+0.5),log10(NgbI[,4]+0.5), col=cpl1(200)[col1], pch=16, cex=0.3)
plot(log10(NgbG[,4]+0.5)/log10(rowSums(Ngba+0.5)),log10(NgbI[,4]+0.5)/log10(rowSums(Ngba+0.5)), col=cpl1(200)[col1], pch=16, cex=0.3)
