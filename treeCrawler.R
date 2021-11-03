treeCrawler=function(L, G, s=1000, Iterations=3, ini=NULL, fl=F, i0=NULL){
  require(igraph)
  g=induced.subgraph(G,L)
  mxc=100
  cl=lapply(1:Iterations, function(i){
    clu=cluster_label_prop(g, initial = ini)
      if (max(clu$membership)==1) return (list(L)) else return(communities(clu))
    })
  l=lengths(cl)
  if (max(l)==1) {
    #print(c(length(L),"single"))
    i=1
    mcl=1
    if (fl==T) {
      iniv=NULL
      mxc=10000
      fl=F
    }
    else iniv=ini    
    while (mcl==1 & i<mxc){
      clu=cluster_label_prop(g, initial = iniv)
      cl=communities(clu)
      mcl=max(clu$membership)
      i=i+1
    }
    if(mcl==1) {
      #print("Failed split")
      fl=T
    }
  }
  else {
    cl=cl[[which.max(l)]] 
    #print(c(max(l),"##",lengths(cl)))
  }
  x=lapply(seq_along(cl),function(i){
    co=cl[[i]]
    print(c("length co ",length(co)))
    if (is.null(ini)) ini1=NULL else {
      x=ini[co]
      ini1=factor(x)
      ini1=as.double(ini1)
      names(ini1)=names(x)
    }
    if (length(co)>=s) {
        y=treeCrawler(co,g,ini=ini1, fl=fl, i0=i0)
        return(y)
      }
      else {
        sink("treemon.txt", append = T)
        print(c(length(co),i0))
        sink()
        return(co)
      }
  })
  return(x)
}