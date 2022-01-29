require(stringdist)
require(parallel)
require(matrixStats)
require(igraph)
require(reshape2)
require(stringi)
require(gplots)
require(ggplot2)
require(qualV)
require(msa)
require(vioplot)
require(Biostrings)
require(seqinr)
require(pbapply)
require(emmeans)
require(Rcapture)
require(eulerr)
require(ggseqlogo)
require(entropy)
require(ppcor)
require(rpart)
require(rpart.plot)
require(ks)

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
plot(Vennall, labels=c("C", "A", "G"), quantities=F, adjust_labels=T, cex=0.5, font=4)

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

load("AdjL_a") # Adjacency list of the combined libraries A and C
load("Ngba") # Neighborhood table of cntr and APLS 
colnames(Ngba)=c("Cntr", "APLS", "All")

Ga=adjL2Glite(AdjL_a) # The combined graph
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

x=Ng1[adj.pNn<0.05,] # Signficant neighborhoods
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
pdf(file="Vulcano_aPLvsContr_u.pdf", height=5,width = 5)
plot(log(((Ng1[,2]+0.5)*pn[1])/((Ng1[,1]+0.5)*pn[2])), -log(pNn), pch=16, col=1+(adj.pNn<0.05)+(log(Ng1[,2]/Ng1[,1])<1&adj.pNn<0.05)*1, cex=-log(pNn)/45+0.01, xlab="log Fold Predominance", ylab="-log(p)", xlim=c(-3,3))
dev.off()

LOR=log10(x[,2]/x[,1])
hist(LOR, breaks=100, xlim=c(-1,1), main="Distribution of log odds of significant neighborhoods")
text(c(-0.75, 0.75), c(20, 20), labels=c("Control", "APLS"))

vulcalls=sort(LOR) # LOR of the A/C proportion in significant neighborhoods, names - central vertex 
save(vulcalls, file="vulcalls")

fln=fla[rownames(Ng1)]
yr=range(-log(pNn))
gr=fln[adj.pNn<0.05]
y=log((Ng1[adj.pNn<0.05,2]*pn[1])/(Ng1[adj.pNn<0.05,1]*pn[2]))
grr=c("Contr.","aPL","Contr.+aPL")
grrs=grr[gr]
vioplot(y[grrs=="Contr."],y[grrs=="aPL"],y[grrs=="Contr.+aPL"], names=c("C","A","C+A"), horizontal = T, ylim=c(-3,3))



# Graph clusterings -------------------------------------------------------

load("Ga")
load("fla")
load("matpep")

dGa=degree(Ga)
vGv=names(V(Ga))
ccGa=transitivity(Ga, type="localundirected")
ccGag=transitivity(Ga, type="global")
ccGa[is.na(ccGa)]=0
save(ccGa, file="ccGa")
save(vGv,file="vGv")
flacl=fla[vGv] # diagnosis label of the vertices of Ga
tfl=table(flacl)


# Louvain clustering ------------------------------------------------------

N=380000
GaLouv=treeCrawLouv(vGv, Ga) #recursive Louvain clustering - yields a tree
maxc=c()
x=GaLouv
repeat{            # finding the cutting level which produces the max N of branches with >3 leaves
  print(length(x))
  x=unlist(x, recursive = F)
  xi=lapply(x, unlist)
  if (length(xi[lengths(xi)==1])>N) break
  n=length(xi[lengths(xi)>3])
  print(n)
  maxc=c(maxc,n)
}
x=GaLouv
j=which.max(maxc)  # tree cut at the max branching level 
for (i in 1:j) x=unlist(x, recursive = F)
x=lapply(x, unlist)
Galouv1_3=unlist(x[lengths(x)<4])
table(flacl[Galouv1_3])/length(Galouv1_3)
xup=unlist(x[lengths(x)>3])
table(flacl[xup])/length(xup)

x=x[lengths(x)>3]  # filtering small clusters
print(length(x))
names(x)=seq_along(x)
x=lapply(x,function(l) {
  names(l)=NULL
  return(l)
})
Galouv=x  # Final clustering of the whole graph - WGcl
save(Galouv,file="Galouv")


# Distribution of Control and APLS mimotopes by cluster -------------------

t0=1:3
names(t0)=t0
ultfl=t(sapply(Galouv, function(l){
  t=t0
  t1=table(flacl[l])
  t[names(t1)]=t1
  return(t)
}))
ultfl2=round(ultfl[,1:2]+ultfl[,3]/2)
pN=cbind(ultfl2[,1], rowSums(ultfl2))
pn=(tfl[1:2]+tfl[3]/2)/sum(tfl) # marginal frequency

cl=makeCluster(4)  # By row proportion/binomial test relative to marginal frequency
clusterExport(cl, "pn")
pNn=pbapply(pN,1,function(l){
  if (l[2]>30) prop.test(l[1], l[2], p=pn[1])$p.value else binom.test(l[1], l[2], p=pn[1])$p.value
}, cl=cl)
stopCluster(cl)

apNn=p.adjust(pNn)
chqult=chisq.test(ultfl, simulate.p.value = T)
rownames(chqult$stdres)=seq_along(chqult$stdres[,1])
y=chqult$stdres[apNn<0.05,]
y=y[order(y[,1]),]
code=c(rep("A",sum(y[,1]<0)), rep("C",sum(y[,1]>0)))

ultimateClCalls=Galouv[as.numeric(rownames(y))]
for (l in ultimateClCalls) print(l[1:10])

CntrClbig=ultimateClCalls[code=="C"] #Clusters in WGcl with underexpression
APLSClbig=ultimateClCalls[code=="A"] #Clusters in WGcl with overexpression
pepCbig=unique(unlist(CntrClbig))
pepAbig=unique(unlist(APLSClbig))
callspep=unlist(ultimateClCalls)

rultfl=log10(ultfl2[,2]*pn[1]/(ultfl2[,1]*pn[2]))

pdf(file="vulcan0_clusters.pdf", width=5, height=5)
plot(rultfl,-log10(pNn), cex=-log10(pNn)/50+0.3, pch=16, col=(apNn<0.05)+1+((rultfl>0)&(apNn<0.05)), xlab="Log ratio APLS to Control", ylab="-log10(p)")
dev.off()


#  Degree Distr and CC by degree plots for Ga -----------------------------

x=cut(dGa,50, labels = F)
y=aggregate(dGa, by=list(x), "length")
colnames(y)=c("Degree","Probability")
z=aggregate(dGa, by=list(x), function(a) 10^(mean(log10(a))))[,2]
y$Degree=z
y$Probability=y$Probability/sum(y$Probability)
xr=range(y$Degree)
yr=range(y$Probability)
plot(y, log="xy", pch=16, xlim=xr, ylim=yr, cex=0.5)
par(new=T)
lmtg=lm(log(y$Probability[6:(length(z))])~y$Degree[6:(length(z))])
x=summary(lmtg)
print(x)
C=exp(x$coefficients[1,1])
print(C)
plot(y$Degree[6:(length(z))],C*exp(x$coefficients[2,1]*y$Degree[6:(length(z))]), 
     log="xy", ty="l",lwd=3,  col=2, xlim=xr, ylim=yr,xlab="", ylab="")
legend("bottomleft", legend = "Exponential curve 0.38*exp(-0.006*k) ", 
       lwd=2, col=2, bty = "n")

snglts=setdiff(names(fla), vGv)
save(snglts, file="snglts")

dccGa=aggregate(ccGa, by=list(dGa), "mean")
cctrnd=mean(dGa)/length(vGv)
pdf(file="CCbyD.pdf", width=5.7,height=6)
plot(dccGa[-(1:2),1],dccGa[-(1:2),2]/cctrnd, log="xy", pch=16, cex=0.5,xlim=c(3,2000),
     xlab="Degree",ylab="Clustering coefficient (fold increase from random graph)")
dev.off()


# Comparison to the significant neighborhood vertices ---------------------

load("vulcalls")
vnm=names(vulcalls)
load("AdjL_a")
egovnm=AdjL_a[vnm]
egovnm=unique(c(unlist(egovnm),names(egovnm)))

Gvulc=induced_subgraph(Ga, egovnm)
vGvu=names(V(Gvulc))
dGvu=degree(Gvulc)
save(Gvulc,file="Gvulc")

GvuLouv=treeCrawLouv(vGvu, Gvulc)

N=0.67*length(vGvu)
maxc=c()
x=GvuLouv
repeat{
  print(length(x))
  x=unlist(x, recursive = F)
  xi=lapply(x, unlist)
  if (length(xi[lengths(xi)==1])>N) break
  n=length(xi[lengths(xi)>3])
  print(n)
  maxc=c(maxc,n)
}
x=GvuLouv
j=which.max(maxc)
for (i in 1:j) x=unlist(x, recursive = F)
x=lapply(x, unlist)
x=x[lengths(x)>3]
print(length(x))
names(x)=seq_along(x)
x=lapply(x,function(l) {
  names(l)=NULL
  return(l)
})
Gvulouv=x
save(Gvulouv,file="Gvulouv")

t0=1:3
names(t0)=t0
ultvufl=t(sapply(Gvulouv, function(l){
  t=t0
  t1=table(fla[l])
  t[names(t1)]=t1
  return(t)
}))
ultvufl2=round(ultvufl[,1:2]+ultvufl[,3]/2)
pNvu=t(apply(ultvufl2,1,function(l){
  c(round(l[1]), sum(l))
}))
pnvu=colSums(ultvufl2)/sum(ultvufl2)
cl=makeCluster(4)
clusterExport(cl, "pnvu",envir = environment())
pNnvu=pbapply(pNvu,1,function(l){
  if (l[2]>30) prop.test(l[1], l[2], p=pnvu[1])$p.value else binom.test(l[1], l[2], p=pnvu[1])$p.value
}, cl=cl)
stopCluster(cl)
apNnvu=p.adjust(pNnvu)
chqultvu=chisq.test(ultvufl, simulate.p.value = T)
rownames(chqultvu$stdres)=seq_along(chqultvu$stdres[,1])
y=chqultvu$stdres[apNnvu<0.05,]
y=y[order(y[,1]),]
codevu=c(rep("A",sum(y[,1]<0)), rep("C",sum(y[,1]>0)))

ultimatevuClCalls=Gvulouv[as.numeric(rownames(y))]
for (l in ultimatevuClCalls) print(l[1:10])

rultvufl=log10(ultvufl2[,2]*pnvu[1]/(ultvufl2[,1]*pnvu[2]))
pdf(file="vulcan0_vuclusters.pdf", width=10, height=10)
plot(rultvufl,-log10(pNnvu), cex=-log10(pNnvu)/30+0.5, pch=16, col=(apNnvu<0.05)+1+((rultvufl>0)&(apNnvu<0.05)), xlab="Log ratio APLS to Control", ylab="Negative logarithm of p")
dev.off()

callsvupep=unique(unlist(ultimatevuClCalls))
allpepcalls1=union(callspep, callsvupep)
fallp1=(allpepcalls1 %in% callspep)*1+(allpepcalls1 %in% callsvupep)*2
names(fallp1)=allpepcalls1
pdf(file="BoxplotDeg1.pdf", width=6, height = 10)
boxplot(log10(dGa[allpepcalls1])~fallp1, notch=T, names=c("Clusters", "Neighborhoods", "Both"), xlab=NULL, ylab="Log degree")
dev.off()

CntrClvu=ultimatevuClCalls[codevu=="C"]
APLSClvu=ultimatevuClCalls[codevu=="A"]
pepCvu=unique(unlist(CntrClvu))
pepAvu=unique(unlist(APLSClvu))

# Combining WGcl and NGcl

names(APLSClbig)=paste(names(APLSClbig),"A", sep="_")
names(APLSClvu)=paste(names(APLSClvu),"Av", sep="_")
names(CntrClbig)=paste(names(CntrClbig),"C", sep="_")
names(CntrClvu)=paste(names(CntrClvu),"Cv", sep="_")
allsigCl=c(CntrClbig,APLSClbig,CntrClvu,APLSClvu)
allsigCl=allsigCl[order(lengths(allsigCl), decreasing = T)]
save(allsigCl, file = "allsigCl")

totallcalls=list(AllCntr=CntrClbig, AllAPLS=APLSClbig, NContr=CntrClvu, NAPLS=APLSClvu)
totallcalls=melt(totallcalls)

nms=names(allsigCl)
dia=strsplit(nms, split="_")
dia=sapply(dia,function(x) x[2])
dia=as.factor(dia)
dian=as.numeric(dia)
D=grepl("A",dia)
szcl=lengths(allsigCl)
szN=lengths(c(CntrClvu, APLSClvu))
dN=c(rep(1,length(CntrClvu)), rep(2,length(APLSClvu)))
szWG=lengths(c(CntrClbig, APLSClbig))
dWG=c(rep(1,length(CntrClbig)), rep(2,length(APLSClbig)))

wilcox.test(szcl~grepl("A",dia)) # APS specific clusters are bigger
wilcox.test(szWG~dWG)
tbl=cbind(table(!names(Galouv) %in% unlist(stri_extract_all(names(CntrClbig), regex="\\d+(?=_)"))), table(!names(Galouv) %in% unlist(stri_extract_all(names(APLSClbig), regex="\\d+(?=_)"))))
prop.test(tbl, alternative = "greater")
wilcox.test(szN~dN)
tbl=cbind(table((names(Gvulouv) %in% unlist(stri_extract_all(names(CntrClvu), regex="\\d+(?=_)")))), table(!names(Galouv) %in% unlist(stri_extract_all(names(APLSClvu), regex="\\d+(?=_)"))))
chisq.test(c(tbl), simulate.p.value = T)

CcallsMx=sapply(CntrClbig, function(Cc){
  sapply(CntrClvu, function(Cv){
    length(intersect(Cc,Cv))/length(unique(c(Cc,Cv)))
  })
})

AcallsMx=sapply(APLSClbig, function(Cc){
  sapply(APLSClvu, function(Cv){
    length(intersect(Cc,Cv))/length(unique(c(Cc,Cv)))
  })
})

# degree stats for significant clusters

degtotc=cbind(totallcalls,dGa[totallcalls[,1]])
tx=table(degtotc$L1)
LGa=length(vGv)
degBS=sapply(tx,function(i){
  sapply(1:1000, function(j){
    median(dGa[sample(LGa,i)])
  })
})
degBS=melt(degBS)
degBS5_95=t(sapply(unique(degBS$Var2), function(x){
  quantile(log10(degBS[degBS$Var2==x,3]),c(0.05,0.95))
}))

pdf(file="BoxplotDegCl.pdf", width=6, height = 10)
boxplot(degtotc[,4]~degtotc[,3], notch=T, xlab=NULL, ylab="Degree", log="y")
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(degBS5_95[,1],each=2)), col=2)
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(degBS5_95[,2],each=2)), col=2)
dev.off()

# eigen centrality stat for significant clusters

eigcGa=eigen_centrality(Ga)
eigtotc=cbind(totallcalls,eigcGa$vector[totallcalls[,1]])
save(eigcGa,file="eigcGa")
eigBS=sapply(tx,function(i){
  sapply(1:1000, function(j){
    median(eigcGa$vector[sample(LGa,i)])
  })
})
eigBS=melt(eigBS)
eigBS5_95=t(sapply(unique(eigBS$Var2), function(x){
  quantile(log10(eigBS[eigBS$Var2==x,3]),c(0.05,0.95))
}))

pdf(file="BoxplotEigCl.pdf", width=6, height = 10)
boxplot(eigtotc[,4]~eigtotc[,3], notch=T, xlab=NULL, ylab="Eigen centrality", log="y")
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(eigBS5_95[,1],each=2)), col=2)
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(eigBS5_95[,2],each=2)), col=2)
dev.off()

# modularity stats for significant cluster

modvuA=Clumodul(Ga,APLSClvu)
modvuC=Clumodul(Ga,CntrClvu)
modClA=Clumodul(Ga,APLSClbig)
modClC=Clumodul(Ga,CntrClbig)
modCl=Clumodul(Ga,Galouv)
modClvu=Clumodul(Ga,Gvulouv)

CLClC=ClmodBS(Ga, Galouv, length(CntrClbig))
f=ecdf(CLClC)
pC=f(modClC)

CLClA=ClmodBS(Ga, Galouv, length(APLSClbig))
f=ecdf(CLClA)
pA=f(modClA)

CLCluC=ClmodBS(Ga, Gvulouv, length(CntrClvu))
f=ecdf(CLCluC)
puC=f(modvuC)

CLCluA=ClmodBS(Ga, Gvulouv, length(APLSClvu))
f=ecdf(CLCluA)
puA=f(modvuA)

paro=par
par(mai=c(1.33, 0.8, 0.8, 0.4)+0.02)
barplot(c(modCl,modClC,modClA,modClvu,modvuC,modvuA), names.arg = c("All", "Control", "APLS", "Neighborhoods","N-Control","N-APLS"), ylab = "Modularity of cluster induced subgraphs", ylim=c(0,1), las=2, yaxt="n")
boxplot(CLClC, notch=T, add=T, at=1.9, ylim=c(0,1), col=0, border=2, cex=0.3)
boxplot(CLClA, notch=T, add=T, at=3.1, ylim=c(0,1), col=0, border=2)
boxplot(CLCluC, notch=T, add=T, at=5.5, ylim=c(0,1), col=0, border=2, cex=0.3)
boxplot(CLCluA, notch=T, add=T, at=6.7, ylim=c(0,1), col=0, border=2)
text(c(1.9,3.1,5.5,6.7), rep(0.95, 4), labels = paste("p=",c(pC,pA,puC,puA), sep=""))
par=paro

CLfla=ClmodBS1(Ga, flacl)
modfla=modularity(Ga,flacl)
f=ecdf(CLfla)
pfla=f(modfla)

hist(CLfla, breaks = 100, xlim = c(-1e-3,6e-3))
lines(c(modfla,modfla),c(0,40), col=2, xlim = c(-1e-3,6e-3))

Dga=diameter(Ga) # too slow

proct=proc.time()
FD=t(sapply(sample(vGv,10000), function(v){
  vi=sapply(1:3, function(i){
    length(ego(Ga,i,v)[[1]])
  })
  c(vi,dGa[v])
}))
print(proc.time()-proct)

x=lm(log10(t(FD[,1:3]))~log10(1:3))$coefficients[2,]
plot(log10(FD[,4]), x)

x=which(dGa<4)
y=distances(Ga, v=x, to=x)
table(y)  # diameter of the graph - 11, 99.3% of the paths are <7 

tdGa=table(dGa)
pdf(file="DegreeDist.pdf",  width = 7, height=7)
# xi=seq(-0.05,3.3,0.05)
# xi=xi[-c(3:8,10,11,13,14,17,20)]
# z=cut(dGa,10^xi, labels=F)
# y1=table(cut(dGa,c(10^xi), labels=F))
# y1=y1/sum(y1)
# x1=aggregate(dGa, by=list(z), function(x) 10^(mean(log10(x))))[,2]        #
# plot(x1,y1, log="xy", ylim=c(1e-5,1e-1), xlab="Degree", ylab="Probability", pch=16)

y1=table(dGa)
x1=as.numeric(names(y1))
plot(x1, y1/sum(y1), log="xy", ylim=c(1e-5,1e-2), xlab="", ylab="", cex=0.5, pch=16)
ytick=c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
axis(side=2, at=ytick, labels = ytick)
dev.off()


#  Mapping global and idiotopes -------------------------------------------

load("NgbI")
load("NgbG")
load("matpep")
NI=rowSums(NgbI[,4:7])
NG=rowSums(NgbG[,4:7])
LNN=length(NI)

cl=makeCluster(3)
ij=cut(seq_along(vGv), 300, labels=F)
clusterExport(cl, c("matpep","vGv","ij"), envir = environment())
clusterEvalQ(cl, require(stringdist))
matmap=pbsapply(1:300, function(i){
  pp=vGv[ij==i]
  sapply(pp,function(p) sum(stringdist(p,matpep, method = "lcs")<5))
}, cl=cl)
stopCluster(cl)
save(matmap, file="matmap")
matmap=unlist(matmap)

NNBS5_95=t(sapply(tx,function(i){
  x=sapply(1:1000, function(j){
    xi=median(NI[sample(LNN,i)])
    xg=median(NG[sample(LNN,i)])
    return(c(xi,xg))
  })
  qi=quantile(log10(x[1,]),c(0.05,0.95))
  qg=quantile(log10(x[2,]),c(0.05,0.95)) 
  return(c(qi,qg))
}))

pdf(file="BoxplotNI.pdf", width=6, height = 10)
boxplot(NI[totallcalls$value]~totallcalls$L1, notch=T, ylab="N of idiotope neighbors", xlab="", log="y")
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(NNBS5_95[,1],each=2)), col=2)
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(NNBS5_95[,2],each=2)), col=2)
dev.off()

pdf(file="BoxplotNG.pdf", width=6, height = 10)
boxplot(NG[totallcalls$value]~totallcalls$L1, notch=T, ylab="N of public mimotope neighbors", xlab="", log="y")
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(NNBS5_95[,3],each=2)), col=2)
lines(c(0,1.5,1.5,2.5,2.5,3.5,3.5,4.5),10^(rep(NNBS5_95[,4],each=2)), col=2)
dev.off()

y=lm(log10(NI[totallcalls$value])~log10(NG[totallcalls$value]))
summary(y)

summary(lm(log10(NG[vGv]+0.5)~log10(eigcGa$vector[vGv])))
ggplot(data.frame(x=log10(eigcGa$vector[vGv]),y=log10(NG[vGv]+0.5)), aes(x=x, y=y) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

summary(lm(log10(NG[vGv]+0.5)~log10(dGa[vGv])))
ggplot(data.frame(x=log10(dGa[vGv]),y=log10(NG[vGv]+0.5)), aes(x=x, y=y) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

summary(lm(log10(NG[vGv]+0.5)~log10(matmap[vGv]+0.5)))
ggplot(data.frame(x=log10(matmap[vGv]+0.5),y=log10(NG[vGv]+0.5)), aes(x=x, y=y) ) +
  geom_bin2d(bins = 100) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()


# comparison of G neighbors in C and A clusters on degree controlling for B
# see LM
# the same for idiotopes
# see LM

x=(lm(log10(dGa)~log10(matmap[vGv]+0.5))$residuals)
plot(log10(dGa),log10((NG[vGv]+0.5))/(x-min(x)+.5), pch=16, cex=0.5, col=rgb(0,0,0,0.1))
plot(log10(dGa),log10((NI[vGv]+0.5))/(x-min(x)+.5), pch=16, cex=0.5, col=rgb(0,0,0,0.1))
plot(log10(dGa),log10(NI[vGv]/NG[vGv]))
mnNIr=sum(NI[vGv])/(sum(NI[vGv])+sum(dGa))

y=cbind(NI[vGv],(NI[vGv]+dGa), dGa)
y=y[order(dGa),]
inor=seq_along(y[,1])[order(dGa)]
sp=500
L=nrow(y)


cl=makeCluster(4)
clusterExport(cl, c("mnNIr", "ij","NIrbyD","dGa","y","sp","L"), envir = environment())
pNNI=pbsapply(seq_along(y[,1]),function(i){
  l=y[i,]
  j1=max(0,i-sp)
  j2=min(L,i+sp)
  j=j1:j2
  pn=sum(y[j,1])/sum(y[j,2])
  if (l[2]>30) prop.test(l[1], l[2], p=pn)$p.value else binom.test(l[1], l[2], p=pn)$p.value
}, cl=cl)
stopCluster(cl)
apNNI=p.adjust(pNNI) # [order(inor)]


plot(log10(dGa),log10(NI[vGv]+0.5)/(NI[vGv]+dGa+.5), pch=16, cex=0.5, col=rgb((apNNI<0.05)*1,0,0,0.1))

Gallcalls=induced_subgraph(Ga, allpepcalls1)
vGall=names(V(Gallcalls))
fClAPLS=1*(allpepcalls1 %in% union(unlist(APLSClbig),unlist(APLSClvu)))+1
names(fClAPLS)=allpepcalls1
Gallcalls=set_vertex_attr(Gallcalls, name = "Dia", value = fla[vGall])
Gallcalls=set_vertex_attr(Gallcalls, name = "nnId", value = log10(NI[vGall]))
Gallcalls=set_vertex_attr(Gallcalls, name = "nnGl", value = log10(NG[vGall]))
Gallcalls=set_vertex_attr(Gallcalls, name = "cluA", value = fClAPLS[vGall])

write.graph(Gallcalls, file="Gallcalls.graphml", format = "graphml")

clsiNI=sapply(allsigCl,function(cl) median(NI[cl]))
clsiNG=sapply(allsigCl,function(cl) median(NG[cl]))



# Peptides from the proteome ----------------------------------------------

aa=AA_ALPHABET[1:20]
proteome=read.fasta("GRCh38_latest_protein.faa", seqtype = "AA")
prots=lapply(proteome, function(p){
  p[p=="U"]="C"
  p[!(p %in% aa)]="A"
  p1=paste(p, collapse="",sep="")
  attributes(p1)=attributes(p)
  return(p1)
})
rm(proteome)
protann=lapply(prots,function(x) attr(x,"Annot"))

cl=makeCluster(4)
clusterExport(cl,"prots")
protsmpep=pbsapply(prots, function(p) {
  i=nchar(p)-6
  sapply(1:i, function(j) substr(p,j,j+6))
}, cl=cl)
stopCluster(cl)

protsmpep=unlist(protsmpep)
protsmpep=unique(protsmpep)  #10.7e6
save(protsmpep, file="proteome7mers")


# PSSM --------------------------------------------------------------------

pssm0=consensusMatrix(vGv, as.prob = T)
x0=array(0,c(20,7))
rownames(x0)=sort(AA_ALPHABET[1:20])

cl=makeCluster(4)
clusterExport(cl,c("x0","pssm0","pscnt"))
clusterEvalQ(cl, require(Biostrings))
allsiPSSM=pblapply(allsigCl,function(l){
  x=consensusMatrix(l)
  N0=length(l)
  x0[rownames(x),]=x
  pscnt(x0,N0,pssm0)
  
}, cl=cl)
stopCluster(cl)

cl=makeCluster(4)
clusterExport(cl,c("allsiPSSM","allsigCl","valign"))
allscopss=pblapply(nms,function(n){
  valign(allsigCl[[n]],allsiPSSM[[n]])  # vectorized alignpssm
}, cl=cl)
stopCluster(cl)
allscopssm=unlist(allscopss)

cl=makeCluster(4)
clusterExport(cl,c("allsiPSSM","protsmpep","valign"))
pepssm=pbsapply(1:52, function(i) {
  valign(sample(protsmpep,1000),allsiPSSM[[i]])
}, cl=cl)
stopCluster(cl)

pepssm=c(pepssm)
pepssmforROC=cbind(c(allscopssm,pepssm),c(rep("1",length(allscopssm)),rep("0",length(pepssm))))
rownames(pepssmforROC)=NULL
pepssmforROC=data.frame(Score=as.numeric(pepssmforROC[,1]), Code=as.numeric(pepssmforROC[,2]))

rocoocpssm=roc(Code~Score, data=pepssmforROC, ci=T)
plot(rocoocpssm, print.auc=T, print.thres=T)
rocoocpssm$auc
rocoocpssm$ci
coords(rocoocpssm, x="best", input="threshold")

ij=seq_along(allsigCl)
cl=makeCluster(4)
clusterExport(cl,c("allsiPSSM","allsigCl","valignpssm","ij"), envir = environment())
clvscl=pbsapply(ij, function(i){
  sapply(ij, function(j){
    valign(allsigCl[[i]],allsiPSSM[[j]])
  })
}, cl=cl)
stopCluster(cl)

plot(log10(szcl), log10(clsiNG))
summary(lm(log10(clsiNI)~log10(szcl)))
summary(lm(log10(clsiNG)~log10(szcl)))



#  Logos of the clusters --------------------------------------------------

allCltit=aggregate(nms, by=list(dia), "list")
pdf(file="clogosA.pdf", width = 8, height = 3.5)
ggplot() + geom_logo(allsigCl[unlist(c(allCltit[1:2,2]))]) + theme_logo() + 
  facet_wrap(~seq_group, ncol=7, scales='free_x')
dev.off()
pdf(file="clogosC.pdf", width = 7, height = 10)
ggplot() + geom_logo(allsigCl[unlist(c(allCltit[3,2]))]) + theme_logo() + 
  facet_wrap(~seq_group, ncol=6, scales='free_x')
dev.off()
pdf(file="clogosCv.pdf", width = 8, height = 12)
ggplot() + geom_logo(allsigCl[unlist(c(allCltit[4,2]))]) + theme_logo() + 
  facet_wrap(~seq_group, ncol=7, scales='free_x')
dev.off()


#  Dominant LCS per cluster check against linear B epis ----------------------

dLCS=lapply(allsigCl, domLCS)
IEDB=read.csv("IEpi.csv")
IEDB$Organism.Name=sapply(seq_along(IEDB$Organism.Name),function(i) {
  x1=IEDB$Organism.Name[i]
  x2=IEDB$Organism.Name.1[i]
  i=nchar(x1)
  x1[i==0]=x2[i==0]
  return(x1)
})
IEDB=IEDB[nchar(IEDB$Organism.Name)>0,]
IEDB=IEDB[,-6]
i=unlist(lapply(strsplit(IEDB$Description,split="") ,function(l){
  all(l %in% AA_ALPHABET[1:20])
}))
IEDB=IEDB[i,] # Linear B cell epitopes 7-25 residues, w/o modifications 
# with known antigen and organism relations

vgrep=Vectorize(grep,vectorize.args = "pattern")

dLCSIEDBhits=lapply(dLCS, function(l){
  i=vgrep(names(l),IEDB$Description )
  lapply(i,function(j) IEDB[j,])
})

dLCShitsA=dLCSIEDBhits[grep("A",dia)]
dLCShitsA=melt(unlist(dLCShitsA, recursive = F))
write.csv(dLCShitsA, file="dLCShitsA.csv")


# Epitope scan by PSSM ----------------------------------------------------

IE7upep=lapply(strsplit(IEDB$Description, split=""),unlist)
L=allsiPSSM[grep("A",dia)]
L=lapply(L,function(M) {
  colnames(M)=1:7
  return(M)
})
cl=makeCluster(4)
clusterExport(cl, "IE7upep", envir = environment())
IEApsscore=pbsapply(L, function(M){
  sapply(IE7upep,function(ep){
    max(sapply(1:(length(ep)-6),function(i) {
      p=ep[i:(i+6)]
      print(ep)
      -log10(prod(M[cbind(p,1:7)]))
    }))
  })    
}, cl=cl)
stopCluster(cl)

jj=IEApsscore<8.3
j=rowSums(jj)
IEAcalls=apply(jj,2,function(i) IEDB$Epitope.ID[i])
IEAcalls=melt(IEAcalls)
IEAcalls=aggregate(IEAcalls$L1, by=list(IEAcalls$value), "list")
IEAcalls=IEAcalls[order(lengths(IEAcalls$x), decreasing = T),]
IEDB[(IEDB$Epitope.ID %in% IEAcalls$Group.1[1:13]),]
write.csv(table(IEDB$Organism.Name[IE7up$Epitope.ID %in% IEAcalls$Group.1]), file="IEDBlinB_APS_PSSMhits.csv")



# Proteome scan ------------------------------------------------------

# Scan all proteome with 52 clusters looking for exact matches of all sequences

cl=makeCluster(4)
clusterExport(cl, c("prots", "vgrep"), envir = environment())
y=allsigCl[order(lengths(allsigCl), decreasing = T)]
prhits=pblapply(y,function(x){
  i=vgrep(x,prots)
  i=i[lengths(i)>0]
  sapply(prots[unlist(i)],function(j) attributes(j)$Annot)
}, cl=cl)
stopCluster(cl)

prh=unlist(prhits, recursive = F)

prhL=unique(unlist(prhits))
length(prhL)/length(prots)
z=lengths(prhits)/(lengths(y))
plot(lengths(y),(z+0.01), col=1+((1:52) %in% grep("A",names(y))), cex=0.75+0.3*((1:52) %in% grep("A",names(y))), pch=16, xlab="Cluster size", ylab="Hits per sequence", log="xy")
l=lm(log10(z+0.01)~log10(lengths(y)))
abline(l)
summary(l)
text(22,1.5, labels="R^2=0.1635,  p=0.00173" )
summary(lm(log10(z+0.01)~log10(clsiNI[order(lengths(allsigCl), decreasing = T)])+log10(szcl)))
(pcor(cbind(log10(z+0.01),log10(clsiNI[order(lengths(allsigCl), decreasing = T)]),log10(szcl)))$estimate)^2

#  Proteome seq hits ------------------------------------

prhP=prots[protann %in% prhL]
prhPann=sapply(prhP,function(x) attributes(x)$Annot)
y=unique(unlist(allsigCl))
yj=cut(seq_along(y),11, labels = F)

cl=makeCluster(4)
clusterExport(cl, c("prhP","prhPann","y","yj"), envir = environment())
clusterEvalQ(cl, require(stringi))

prhitsp=pblapply(1:11,function(i){
  x=y[yj==i]
  a=sapply(x, function(p){
    i=stri_locate_all(prhP, fixed=p)
    names(i)=prhPann
    i[sapply(i,function(l) !is.na(l[[1]]))]
  })
  a[lengths(a)>0]
}, cl=cl)
stopCluster(cl)
prhitsp=unlist(prhitsp, recursive = F)
prhitspp=sapply(prhitsp,function(l) unique(names(l)))
prhitpp=unlist(sapply(seq_along(prhitsp), function(l) names(prhitsp[[l]])))

# 1072 unique mimotope sequences are found in the proteome
prhpp=unique(names(prhitspp))
NIprpp=names(NI) %in% prhpp

pdf(file="BoxplotNIprpp.pdf", width=3, height = 10)
boxplot((NI+1)~NIprpp, notch=T, log="y", xlab="Self peptidome", ylab="N  of idiotype nearest neighbors")
dev.off()

NGprpp=names(NG) %in% prhpp

pdf(file="BoxplotNGprpp.pdf", width=3, height = 10)
boxplot((NG+1)~NGprpp, notch=T, log="y", xlab="Self peptidome", ylab="N  of public mimotopes nearest neighbors")
dev.off()

summary(lm(log10(NI[vGv]+0.5)~(NIprpp[names(NI) %in% vGv]*1+1)+log10(dGa)))
summary(lm(log10(NG[vGv]+0.5)~(NGprpp[names(NG) %in% vGv]*1+1)+log10(dGa)))
summary(lm((log10(dGa)~NIprpp[names(NG) %in% vGv]*1+1)))
boxplot(log10(dGa)~NIprpp[names(NG) %in% vGv]*1+1, notch=T)


prhpprot=lapply(prhitspp, function(n) {
  j=grepl("isoform",n)
  x1=stri_extract_all(n[j], regex="(?<=>.{1,15} ).+(?=iso)")
  x2=stri_extract_all(n[!j], regex="(?<=>.{1,15} ).+(?=\\[Homo)")
  x=c(x1,x2)
  x=x[!is.na(x)]
  unique(x)
})

prhpprot=unlist(prhpprot) # 1367 unique proteins (isoforms are pooled)
write.csv(prhpprot,file="proteinswithhits.csv")



# Entropy -----------------------------------------------------------------

mimsaa=t(sapply(strsplit(vGv,split = ""), unlist))
callsaa=lapply(allsigCl, function(l) t(sapply(strsplit(l,split = ""), unlist)))
protsaa=t(sapply(strsplit(protsmpep,split = ""), unlist))

Ecalls=lapply(callsaa,function(m){
  proct=proc.time()
  print(nrow(m))
  a=apply(m,1,function(l) entropy(table(l), method = "MM"))
  print(proc.time()-proct)
  return(a)
})

proct=proc.time()
Emims=apply(mimsaa,1,function(l) entropy(table(l), method = "MM"))
print(proc.time()-proct)

proct=proc.time()
j=cut(seq_along(protsmpep), 120, labels=F)
cl=makeCluster(4)
clusterExport(cl, c("protsaa","j"))
clusterEvalQ(cl, require(entropy))
Eprots=pbsapply(1:120, function(i){
  apply(protsaa[j==i,],1,function(l) entropy(table(l), method = "MM"))
},cl=cl)
stopCluster(cl)
print(proc.time()-proct)
Eprots=unlist(Eprots)

ecalls=melt(Ecalls)
allE=list(Mimotopes=Emims, Calls=Ecalls, Proteome=Eprots)
allE=melt(allE)
allE$value=cut(allE$value, 3, breaks=c(-0.001,1.57,2,2.2,2.5), labels = F)
tbl=table(allE$value,allE$L1)
tbls=sweep(tbl,2,colSums(tbl),"/")
chsq=chisq.test(tbl, simulate.p.value = )
chsq$stdres
cols=c(rgb(0,0,0,1),rgb(0.5,0.5,0.5,1),rgb(0.75,0.75,0.75,1),rgb(0.95,0.95,0.95,1))
barplot(chsq$stdres, beside=T, col=cols, main="Distribution of peptides by sequence entropy", ylab = "Standardized residuals")
legend("topleft", c("very low","low","high","very high"),fill=cols, bty = "n")
barplot(tbls, beside=T, col=cols, main="Distribution of peptides by sequence entropy", ylab = "Proportion of sequences", log="y")
legend("topleft", c("very low","low","high","very high"),fill=cols, bty = "n")

allcls=allE[!is.na(allE$L2),]
dia=seq_along(allcls[,1]) %in% grep("A", allcls$L2)
allcls=cbind(allcls, Diagnosis=rep("Control",nrow(allcls)))
allcls$Diagnosis[dia]="AFLS"              
chsqcls=chisq.test(table(allcls[,c(1,4)]))
barplot(chsqcls$stdres, beside=T, col=cols, main="Distribution of calls by sequence entropy", ylab = "Standardized residuals")
legend("topleft", c("very low","low","high","very high"),fill=cols, bty = "n")
tbl2=table(allcls[,c(1,4)])
tbl2=sweep(tbl2,2,colSums(tbl2),"/")
barplot(tbl2, beside=T, col=cols, main="Distribution of calls by sequence entropy", ylab = "Proportion of sequences")
legend("topleft", c("very low","low","high","very high"),fill=cols, bty = "n")

entropyexmpls=data.frame(Seq=unlist(allsigCl),Entropy=unlist(Ecalls), stringsAsFactors = F)
entropyexmpls=entropyexmpls[sample(nrow(entropyexmpls),1000),]
entropyexmpls=entropyexmpls[order(entropyexmpls[,2]),]
write.csv(entropyexmpls, file="entropyexmpls.csv")


chsqclust=chisq.test(table(allcls[,c(1,2)]), simulate.p.value = T)
dia=grepl("A", colnames(chsqclust$stdres))
j=colSums(chsqclust$stdres[1:2,])-colSums(chsqclust$stdres[3:4,])
pdf(file="EntropyClus.pdf",width=12, height=30) #, 
par(mar=c(5,10,3,1))
barplot(chsqclust$stdres[,order(j,decreasing = T)], beside=T, cex.axis=2,horiz=T, col=cols, main="Distribution of calls clusters by sequence entropy", xlab = "Standardized residuals", las=2)
legend("topleft", c("very low","low","high","very high"),fill=cols, bty = "n")
par(mar=c(5,2,2,1))
dev.off()


coff=c(-0.001,1.57,2,2.2,2.5)
lab=c("very low","low","high","very high")
#lab=c("6+1","5+2","4+3","5+2x1","4+2+1","2x3+1","3+2x2","4+3x1","3+2+2x1","3x2+1","3+4x1","2x2+3x1","2+5x1","7x1")
EF=cut(Emims,coff,labels = lab)
boxplot(log10(eigcGa$vector)~EF, xlab="Entropy", ylab = "Log10 Eigenvector Centrality", notch=T)
boxplot(log10(dGa)~EF,xlab="Entropy", ylab = "Log10 degree", notch=T)
boxplot(log10(rowSums(NgbG[vGv,4:7]))~EF,xlab="Entropy", ylab = "Log10 N of global repertoire neighbors", notch=T)
boxplot(log10(rowSums(NgbI[vGv,4:7]))~EF,xlab="Entropy", ylab = "Log10 N of idiotope neighbors", notch=T)
boxplot(log10((NgbI[vGv,4:7]+0.3)/(rowSums(NgbI[vGv,c(1:3,5:7)])+0.3))~EF, xlab="Entropy", ylab = "Log10 nearest neighbor ratio of idiotopes to mimtopes", notch=T)
boxplot(log10((rowSums(NgbI[vGv,4:7])+0.3)/(rowSums(NgbG[vGv,4:7])+0.3))~EF, xlab="Entropy", ylab = "Log10 nearest neighbor ratio of idiotopes to global mimtopes", notch=T)
boxplot(log10((NgbG[vGv,4:7]+0.3)/(rowSums(NgbG[vGv,c(1:3,5:7)])+0.3))~EF, xlab="Entropy", ylab = "Log10 nearest neighbor ratio of global to groups sepcific mimtopes", notch=T)
boxplot(log10(ccGa+0.3)~EF, xlab="Entropy", ylab = "Log10 clustering coefficient", notch=T)

GLvpep=unlist(Galouv)
GvLvpep=unlist(Gvulouv)
nonclpep=vGv[!vGv %in% union(GLvpep,GvLvpep)]

names(EF)=vGv
EFs=as.character(EF)
names(EFs)=vGv
EFs[Emims<1]="Low E"
EFs=as.factor(EFs)
tblGL=sapply(list(GLvpep,GvLvpep,nonclpep),function(L) table(EF[L]))
colnames(tblGL)=c("Large clusters","Neighborhood clusters","Non-clustered")
chsq.GL=chisq.test(tblGL, simulate.p.value = T)
barplot(t(chsq.GL$stdres), beside = T, ylab = "Standardized residuals")
legend("bottom", colnames(tblGL), fill=cols, bty="n")
tblGLs=sweep(tblGL,2,colSums(tblGL),"/")
barplot(tblGLs, beside = T, ylab = "Proportion of sequences")
legend("topleft", c("very low","low","high","very high"), fill=cols, bty="n")

names(Emims)=vGv
EmG=Emims[GLvpep]
mGl=melt(Galouv)
x=mGl$L1
names(x)=mGl$value
mGl=x
GLEm=aggregate(EmG, by=list(mGl[names(EmG)]), "mean")
hist(GLEm$x)        

EmGv=Emims[GvLvpep]
mGvl=melt(Gvulouv)
x=mGvl$L1
names(x)=mGvl$value
mGvl=x
GvLEm=aggregate(EmGv, by=list(mGvl[names(EmGv)]), "mean")
cGvLEm=cut(GvLEm$x,3, breaks=c(-0.001,1.57,2,2.2,2.5), labels = F)
plot(table(cGvLEm))        

EFpr=cut(Eprots,coff,labels = lab)
protsmpli=sapply(lab, function(l){
  i=which(EFpr==l)
  sample(i,0.05*length(i))
})
protsmpli=unlist(protsmpli)
protpepsmpl=protsmpep[protsmpli]
save(protpepsmpl, file="protpepsmpl")


#
#  Check for specificity of low complexity patterns in significant --------
#

sx=sample(Galouv[!names(Galouv) %in% nms], 50)
sx=c(sx,allsigCl)

smplPSSM=lapply(sx,function(l){
  x=consensusMatrix(l)
  N0=length(l)
  x0[rownames(x),]=x
  pscnt(x0,N0,pssm0)
})
pdf(file="clogossmp1.pdf", width = 8, height = 50)
ggplot() + geom_logo(sx) + theme_logo() + 
  facet_wrap(~seq_group, ncol=7, scales='free_x')
dev.off()

load("protpepsmpl") 
rm(protsmpep)

cl=makeCluster(4)
clusterExport(cl, c("protpepsmpl","valign"),envir = environment())
sxalpssm=pbsapply(smplPSSM,function(m) {
  valign(protpepsmpl,m)
}, cl=cl)
stopCluster(cl)

save(sxalpssm, file="sxalpssm")

x=strsplit(protpepsmpl, split="")
x=lapply(x,table)

Eprsmp=sapply(x,entropy, method="MM")
names(Eprsmp)=protpepsmpl
save(Eprsmp,file="Eprsmp")
#coff=c(sort(unique(Eprsmp))-0.00001,3)
#lab=c("6+1","5+2","4+3","5+2x1","4+2+1","2x3+1","3+2x2","4+3x1","3+2+2x1","3x2+1","3+4x1","2x2+3x1","2+5x1","7x1")
EFpr=cut(Eprsmp,breaks=c(-0.001,1.57,2,2.2,2.5), labels = F)

tnxalpssm=apply(sxalpssm,2,function(l){
  aggregate(l,by=list(Eprsmp),function(x) sum(x<8.3))
})

tnxalmeans=sapply(tnxalpssm,function(l) l[,2])
fe=cut(tnxalpssm[[1]][,1],breaks=c(-0.001,1.57,2,2.2,2.5), labels = F)
tnxalmeans=t(tnxalmeans)
tnxalmeans=melt(tnxalmeans)
j=rep("Rnd",1428)
j[grep("_",tnxalmeans$Var1)]="Calls"
tnxalmeans[,1]=j
colnames(tnxalmeans)=c("Set","Entropy","N")
tnxalmeans$Entropy=cut(tnxalmeans$Entropy, breaks=c(0.5,1.5,2.5,3.5,4.5), labels = F)
par0=par
par(mai=c(2,1,1,0.5))
boxplot((tnxalmeans$N+0.3)~tnxalmeans$Entropy+tnxalmeans$Set, notch=T, ylab="N positive hits", names=c("A/C very low E","A/C low E","A/C high E","A/C very high E","Rnd very low E","Rnd low E","Rnd high E","Rnd very high E"), xlab=NULL, las=2, log="y")
par=par0


#  LM ----------------------------------------------------------------

load("Emims")
pepCvuneat=pepCvu[!(pepCvu %in% c(pepCbig,pepAbig))]
pepAvuneat=pepAvu[!(pepAvu %in% c(pepCbig,pepAbig))]
C=rep("N", length(vGv))
C[vGv %in% c(pepCbig,pepAbig,pepCvuneat,pepAvuneat)]="C"
C=as.factor(C)
C=relevel(C, "N")
Cdia=rep("N", length(vGv))
Cdia[vGv %in% c(pepCbig,pepCvuneat)]="C"
Cdia[vGv %in% c(pepAbig,pepAvuneat)]="A"
Cdia=as.factor(Cdia)
Cdia=relevel(Cdia, "N")

Cmod=rep("N", length(vGv))
Cmod[vGv %in% c(pepCvuneat,pepAvuneat)]="V"
Cmod[vGv %in% c(pepCbig,pepAbig)]="B"
Cmod=as.factor(Cmod)
Cmod=relevel(Cmod, ref="N")

P=rep("N", length(vGv))
P[vGv %in% prhpp]="P"
P=as.factor(P)
P=relevel(P, ref="N")

bigtbl=data.frame(Calls=C, Dia=Cdia, Cmod=Cmod, P=P, D=log10(dGa+0.5), G=log10(NG[vGv]+0.5), Id=log10(NI[vGv]+0.5), Bg=log10(matmap[vGv]+0.5), E=Emims[vGv])

# linearity
ggplot(bigtbl[sample(1:nrow(bigtbl), 10000),], aes(Bg, D)) +                  
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")

ggplot(bigtbl[sample(1:nrow(bigtbl), 10000),], aes(Bg, Id)) +                 
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")

ggplot(bigtbl[sample(1:nrow(bigtbl), 10000),], aes(D, Id)) +                  
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")

sink(file="pcor.txt")
pcor(bigtbl[,c(5:8)])
sink()

lm=list()
lm[[1]]=summary(lm(data=bigtbl, D~Bg))
lm[[2]]=summary(lm(data=bigtbl, G~Bg))
lm[[3]]=summary(lm(data=bigtbl, Id~Bg))

DBg=lm(data=bigtbl, D~Bg)$residuals
GBg=lm(data=bigtbl, G~Bg)$residuals

smalltbl=bigtbl
smalltbl$Dia=factor(smalltbl$Dia, levels = c("N","C","A"))
smalltbl$G=GBg
smalltbl$D=DBg

GDBg=lm(data=smalltbl, G~D)
lm[[4]]=summary(GDBg)
GDBg=GDBg$residuals
smalltbl=cbind(smalltbl, Gd=GDBg)

IdBg=lm[[3]]$residuals
xt3=cbind(bigtbl,IdBg)
ggplot(xt3[sample(1:nrow(xt3), 10000),], aes(D, IdBg)) +                  
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")
xt3=cbind(smalltbl,IdBg)
ggplot(xt3[sample(1:nrow(xt3), 10000),], aes(D, IdBg)) +                  
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")
rm(xt3)

smalltbl$Id=IdBg
ggplot(smalltbl[sample(1:nrow(smalltbl), 10000),], aes(G, Id)) +                  
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  geom_smooth(method='lm', formula= y~x, color="red")+
  geom_smooth(method='loess', formula= y~x, color="green")

lm[[5]]=summary(glm(data=smalltbl,Dia~(Id+G)^2, family="binomial", subset = smalltbl$Dia!="A"))
lm[[6]]=summary(glm(data=smalltbl,Dia~(Id+G)^2, family="binomial", subset = smalltbl$Dia!="C"))
lm[[7]]=summary(glm(data=smalltbl,Dia~(Id+G)^2, family="binomial", subset = smalltbl$Dia!="N"))
lm[[8]]=summary(glm(data=smalltbl,Calls~(Id+G)^2, family="binomial"))
lm[[9]]=summary(glm(data=smalltbl,P~(Id+G)^2, family="binomial", subset=Dia=="A"))
lm[[10]]=summary(glm(data=smalltbl,P~(Id+G)^2, family="binomial", subset=Dia=="C"))

sink(file="ACNCllsglm.txt")
for (l in lm) print(l)
sink()

mse=aggregate(smalltbl[,6:7], by=list(smalltbl$Dia), function(x) c(mean=mean(x), se=10*sd(x)/sqrt(length(x))))
mse=as.data.frame(cbind(mse[,2], mse[,3]))
colnames(mse)=c("Gm", "Gse","Idm","Idse")
mse=cbind(mse,Dia=levels(smalltbl$Dia))

GIdmapN=kde2d(smalltbl$G[smalltbl$Dia=="N"], smalltbl$Id[smalltbl$Dia=="N"], n=300,lims=c(-0.5,0.5,-1,1))

GIdmapN=kde(smalltbl[smalltbl$Dia=="N",6:7])
GIdmapC=kde(smalltbl[smalltbl$Dia=="C",6:7])
GIdmapA=kde(smalltbl[smalltbl$Dia=="A",6:7])


pdf(file="smalltblIdxG.pdf", width=7, height=7)
plot(GIdmapN,cont=c(10,90), xlim=c(-0.5,0.5),ylim=c(-1,1.5), xlab="G", ylab="Id")
par(new=T)
plot(GIdmapC, cont=c(10,90), xlim=c(-0.5,0.5),ylim=c(-1,1.5), col=rgb(0,1,0,0.8), xlab="", ylab="")
par(new=T)
plot(GIdmapA, cont=c(10,90), xlim=c(-0.5,0.5),ylim=c(-1,1.5), col=rgb(1,0,0,0.8), xlab="", ylab="")
par(new=T)
plot(mse[,c(1,3)],pch=16, col=c(1,2,3), xlim=c(-0.5,0.5),ylim=c(-1,1.5), xlab="", ylab="")
arrows(x0=c(mse$Gm-mse$Gse,mse$Gm),y0=c(mse$Idm,mse$Idm-mse$Idse),x1=c(mse$Gm+mse$Gse,mse$Gm),y1=c(mse$Idm,mse$Idm+mse$Idse), code=3, angle=90, length=0.02, col=c(1,2,3,1,2,3))
abline(lm(data=smalltbl, Id~G), col=4)
#par(new=T)
#plot(smalltbl$G[smalltbl$P=="P"],smalltbl$Id[smalltbl$P=="P"],pch=16, col=(smalltbl$Dia[smalltbl$P=="P"]=="A")*2+(smalltbl$Dia[smalltbl$P=="P"]=="C")*3, cex=0.5, xlim=c(-0.5,0.5),ylim=c(-1,1.5), xlab="", ylab="")
legend("topleft", c("Unchanged", "Enriched in APS", "Lost in APS"), fill=1:3, bty="n")
dev.off()


# Assortativities ---------------------------------------------------------

asDia=assortativity_nominal(Ga,fla[vGv], directed = F)
asD=assortativity(Ga,dGa[vGv], directed = F)
asI=assortativity(Ga,NI[vGv], directed = F)
asG=assortativity(Ga,NG[vGv], directed = F)
asE=assortativity(Ga,Emims[vGv])

ij=cut(rank(dGa), 20, labels = F)
asdg=t(sapply(1:20, function(i){
  G=induced.subgraph(Ga, names(dGa[ij==i]))
  dg=assortativity(G, dGa[names(V(G))])
  ng=assortativity(G, NG[names(V(G))])
  ni=assortativity(G, NI[names(V(G))])
  return(c(Degree=dg,NNpubl=ng,NNIdio=ni))
}))
asall=c(asD,asG, asI)
dx=sapply(sort(unique(ij)),function(i) exp(mean(log(dGa[ij==i]))))
for (i in 1:3){
  plot(dx,asdg[,i], ylim=c(0,1), xlim=c(15,1010),ty="b",col=rgb(i==2,i==3,0,1), pch=16, cex=0.75, xlab="Degree", ylab="Assortativity", log="x")
  par(new=T)
  plot(c(15,1010),c(asall[i],asall[i]), ty="l", lwd=2, col=i, cex=0.75,ylim=c(0,1), xlim=c(15,1010), xlab="", ylab="", log="x")
  par(new=T)
}
plot(c(15,1010),c(0.3,0.3), ty="l", lwd=1, col=4,cex=0.75, ylim=c(0,1), xlim=c(15,1010), xlab="", ylab="", log="x")
par(new=F)
legend("topleft", legend=c("Degree", "NN public", "NN idiotope"), fill=1:3, bty="n")

x1=log10(dGa[vGv]+0.5)
x2=log10(matmap[vGv]+0.5)
y=log10(NG[vGv]+0.5)
summary(lm(x1~x2))
dxGcntrmat=lm(y~x1*x2)
pcor(cbind(x1,x2,y))
(pcor(cbind(x1,x2,y))$estimate)^2

assortativity_nominal(Gvulc,fla[vGvu], directed = F)
assortativity(Gvulc,dGvu[vGvu], directed = F)
assortativity(Gvulc,NI[vGvu], directed = F)
assortativity(Gvulc,NG[vGvu], directed = F)


# J regions most like calls mimotopes --------------------------------------

load("IgJt7")
length(IgJt7)
Acalls=totallcalls$value[grep("APLS",totallcalls$L1)]
IdAdjL=adjL(L1=Acalls,L2=IgJt7)
x=allscopssm[names(IdAdjL)]
plot(x, log10(lengths(IdAdjL)[names(x)]), cex=1.5, pch=16, col=rgb(0,0,0,0.3))
x1=x[x<5&lengths(IdAdjL)[names(x)]>3000]
x2=x[lengths(IdAdjL)[names(x)]>6500]
x=c(x1,x2)
x=IdAdjL[names(x)]
rm(IgJt7)
load("IgJT")
rm(lIgJt)
x=unlist(x)
x=unique(x)
bestAPLcallsIdNN=x
save(bestAPLcallsIdNN, file="bestAPLcallsIdNN")
load("bestAPLcallsIdNN")

IgJtrim=unique(IgJtrim)
actualJs=sapply(seq_along(bestAPLcallsIdNN),function(i){
  x=bestAPLcallsIdNN[[i]]
  IgJtrim[grepl(x,IgJtrim)]
})
actualJs=table(unlist(actualJs))
vgrepl=Vectorize(grepl, vectorize.arg="pattern")
x=bestAPLcallsIdNN[1:100]

j=seq_along(Acalls)
ij=cut(j, 12, labels = F)
cl=makeCluster(4)
clusterExport(cl, c("Acalls","IgJtrim","ij","vgrepl"), envir = environment())
AIdent=pbsapply(1:12, function(i){
  x=Acalls[ij==i]
  apply(vgrepl(x, IgJtrim, fixed=T),2,function(z) IgJtrim[z])
}, cl=cl)
stopCluster(cl)

AIdent=unlist(AIdent, recursive = F)
AIdent=AIdent[lengths(AIdent)>0]
AIdent=unlist(AIdent)
AIdent=unique(AIdent)


load("IgJt7")
Ccalls=totallcalls$value[grep("ntr",totallcalls$L1)]
IdCAdjL=adjL(L1=Ccalls,L2=IgJt7)
x=allscopssm[names(IdCAdjL)]
plot(x, log10(lengths(IdCAdjL)[names(x)]), cex=1.5, pch=16, col=rgb(0,0,0,0.3))
x1=x[x<3.38&lengths(IdCAdjL)[names(x)]>3000]
x2=x[lengths(IdCAdjL)[names(x)]>70000]
x=c(x1,x2)
x=IdCAdjL[names(x)]
rm(IgJt7)
load("IgJT")
rm(lIgJt)
x=unlist(x)
x=unique(x)
bestCntrcallsIdNN=x
save(bestCntrcallsIdNN, file="bestCntrcallsIdNN")
load("bestAPLcallsIdNN")

j=seq_along(Ccalls)
ij=cut(j, 40, labels = F)
cl=makeCluster(3)
clusterExport(cl, c("Ccalls","IgJtrim","ij","vgrepl"), envir = environment())
CIdent=pbsapply(1:40, function(i){
  x=Ccalls[ij==i]
  x=apply(vgrepl(x, IgJtrim, fixed=T),2,function(z) IgJtrim[z])
  x=unique(unlist(x))
  save(x,file=paste("Cident_",i,sep="", collapse=""))
  return(NULL)
}, cl=cl)
stopCluster(cl)


CIdent=sapply(1:40, function(i){
  load(paste("Cident_",i,sep="", collapse=""))
  return(x)
})
CIdent=unique(unlist(CIdent))
save(CIdent,file = "CIdent")
y1=table(nchar(CIdent))
y2=table(nchar(AIdent))
x1=as.numeric(names(y1))
x2=as.numeric(names(y2))
y3=table(nchar(IgJtrim))
x3=as.numeric(names(y3))
plot(x1,y1, xlim=c(0,40), col=rgb(1,0,0,0.5), ty="b", xlab="CDR3 Length")
par(new=T)
plot(x2,y2, xlim=c(0,40), col=rgb(0,0,0,0.5), ty="b", xlab="")
par(new=T)
plot(x3,y3, xlim=c(0,40), col=rgb(0,1,0,0.5), ty="b", xlab="")

wilcox.test(nchar(IgJtrim), nchar(CIdent))
wilcox.test(nchar(IgJtrim), nchar(AIdent))
wilcox.test(nchar(AIdent), nchar(CIdent))
mean(nchar(IgJtrim))
mean(nchar(CIdent))
mean(nchar(AIdent))
mean(nchar(c(AIdent,CIdent)))

aaCId=AAStringSet(CIdent)
aafCId=alphabetFrequency(aaCId)[,1:20]
aafCId=colSums(aafCId)/sum(aafCId)

aaAId=AAStringSet(AIdent)
aafAId=alphabetFrequency(aaAId)[,1:20]
aafAId=colSums(aafAId)/sum(aafAId)


chsqACId=chisq.test(cbind(colSums(alphabetFrequency(aaCId)[,1:20]),colSums(alphabetFrequency(aaAId)[,1:20])), simulate.p.value =T)
chsqACId$stdres
cols=c(1,2,1,2,1,2,1,2,8,0,8,0,8,0,8,0,8,0,1,2,8,0,8,0,8,0,8,0,1,2,1,2,1,2,8,0, 1,2,1,2)
barplot(rbind(aafCId,aafAId), beside = T, ylab="Frequency", col=cols)
legend("topleft", legend = c("Significantly different - Control","Significantly different - APLA","Control","APLA"), fill=c(1,2,8,0), bty="n")

