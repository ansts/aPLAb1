require(XML)
require(stringi)
require(factoextra)
require(NbClust)
require(Rtsne)
require(ggseqlogo)
require(png)
require(imager)
require(Biostrings)
require(ggplot2)
mim7lindb=xmlToList(xmlParse(file="D:\\Documents\\BG\\FNI2016\\work\\PavlinaPap\\mimodb\\AdvanceSearch_By_monoclonal antibody(TargetType)+7(LibrarySeqLength)+Linear(LibraryTopology).xml"), simplify = T)
m7lpep=lapply(mim7lindb[[1]], function(x){
  y=unlist(stri_extract_all_regex(x$Peptides,"[A-Z]\\w+"))
  y=y[sapply(y,nchar)==7]
  return(unlist(y))
})
m7lpep=m7lpep[lengths(m7lpep)>0]
names(m7lpep)=NULL

m7lpbad=as.character(read.table("m7lbad.txt")[,1])
m7lpeps=lapply(m7lpep,function(l){
  l=l[!l %in% m7lpbad]
})

Xs=sapply(m7lpeps,function(s){
  sapply(s,function(p){
    x=stri_locate(p,fixed = "X")
    print(x)
    if (all(is.na(x))) x=0 else x=x[1]
    return(x)
  })
})
m7lpeps=lapply(Xs,function(l){names(l[l==0])})

m7lpeps=m7lpeps[lengths(m7lpeps)>2]

m7lpeps=m7lpeps[-48]
m7lpeps[[49]]=union(m7lpeps[[48]],m7lpeps[[49]])
m7lpeps=m7lpeps[-48]

m7lpeps[[14]]=union(m7lpeps[[14]],m7lpeps[[28]])
m7lpeps=m7lpeps[-28]

write(unlist(m7lpeps), "m7lp.txt")

L7=melt(m7lpepsc)

pdup=L7[duplicated(L7[,1]),1]
L7c=L7[!(L7[,1] %in% pdup),]
tn=table(L7c[,2])
tnc=as.double(names(tn[tn>1]))
L7c=L7c[L7c[,2] %in% tnc,]