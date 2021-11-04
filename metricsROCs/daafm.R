# A string metric based comparison between AA sequences taking into consideration 
# evolutionary conservation. The groups of homologous AA represent clusters of AA
# with all positive scores inbetween in BLOSUM62.

daafm=function(lp1,lp2,mth="lcs", q=2){
  require(Biostrings)
  require(parallel)
  require(stringdist)
  
  aa=AA_ALPHABET[1:20]
  aafmx=matrix("X",20,20)
  rownames(aafmx)=aa
  colnames(aafmx)=aa
  # Grouping the aa in groups which have only positive scores in BLOSUM62 within the group
  aafm=list(M=c("M","I","L","V"), "F"=c("F","W","Y"), Y=c("Y","H"), D=c("D","E"), N=c("D","N"), H=c("H","N"), S=c("S","N"), E=c("E","Q"), K=c("E","K"), R=c("R", "K","Q"), A=c("S","A"), T=c("T","S"), G="G", C="C",P="P")
  for (i in names(aafm)){
      aafmx[aafm[i][[1]],aafm[i][[1]]]=i
  }
  
  dp=sapply(lp1, function(p1){
    p1=unlist(strsplit(p1, split=""))
    cp1=tolower(p1)
    dp0=sapply(lp2, function(p2){
      p2=unlist(strsplit(p2, split=""))
      cp2=tolower(p2)
      m=aafmx[p1,p2]
      n1=length(p1)
      n2=length(p2)
      #if (sum(m!="X")>3) {
        ccp1=sapply(1:n1, function(i){ 
          nx=m[i,]!="X"
          if (sum(nx)>0) {
            return(unique(m[i,nx]))
            }
          else {
            return(cp1[i])
          }
        })
        ccp2=sapply(1:n2, function(i){ 
          nx=m[,i]!="X"
          if (sum(nx)>0) {
            return(unique(m[nx,i]))
          }
          else {
            return(cp2[i])
          }
        })
        ccpL1=c()
        i=1
      for (i1 in 1:length(ccp1[[1]])){
        for (i2 in 1:length(ccp1[[2]])){
          for (i3 in 1:length(ccp1[[3]])){
            for (i4 in 1:length(ccp1[[4]])){
              for (i5 in 1:length(ccp1[[5]])){
                for (i6 in 1:length(ccp1[[6]])){
                  for (i7 in 1:length(ccp1[[7]])){
                    ccpL1[[i]]=c(ccp1[[1]][i1],ccp1[[2]][i2],ccp1[[3]][i3],ccp1[[4]][i4],ccp1[[5]][i5],ccp1[[6]][i6],ccp1[[7]][i7])
                    i=i+1
                  }
                }
              }
            }
          }  
        }
      }
      ccpL1=sapply(ccpL1, function(p){paste(p, collapse="", sep = "")})  
      ccpL2=c()
      i=1
      for (i1 in 1:length(ccp2[[1]])){
        for (i2 in 1:length(ccp2[[2]])){
          for (i3 in 1:length(ccp2[[3]])){
            for (i4 in 1:length(ccp2[[4]])){
              for (i5 in 1:length(ccp2[[5]])){
                for (i6 in 1:length(ccp2[[6]])){
                  for (i7 in 1:length(ccp2[[7]])){
                    ccpL2[[i]]=c(ccp2[[1]][i1],ccp2[[2]][i2],ccp2[[3]][i3],ccp2[[4]][i4],ccp2[[5]][i5],ccp2[[6]][i6],ccp2[[7]][i7])
                    i=i+1
                  }
                }
              }
            }
          }  
        }
      }
      ccpL2=sapply(ccpL2, function(p){paste(p, collapse="", sep = "")})  
      d=sapply(ccpL1, function(p1){
        d0=sapply(ccpL2, function(p2){
          stringdist(p1,p2, method=mth, q=q)
        })
      })
      d=min(d)
      return(d)
    })
    return(dp0)
  })
  
  return(dp)
  
}