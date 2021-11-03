treeCr=function(L, G, s=1000, Iterations=3, ini=NULL){
  require(igraph)
  g=induced.subgraph(G,L)
  if (is.null(ini)) mxc=100 else mxc=5
  cl=lapply(1:Iterations, function(i){
    clu=cluster_label_prop(g, initial = ini)
    if (max(clu$membership)==1) return (list(L)) else return(communities(clu))
  })
  l=lengths(cl)
  #if (max(l)==1) {
    #print(c(length(L),"single"))
    #i=1
    #mcl=1
    #while (mcl==1 & i<mxc) {
    #  clu=cluster_label_prop(g)
    #  cl=communities(clu)
    #  mcl=max(clu$membership)
    #  i=i+1
    #}
    #print(c("mcl ",mcl,"n=",i-1))
  #}
  #else {
    cl=cl[[which.max(l)]]
    #print(c(max(l),"##",lengths(cl)))
  #}
  x=lapply(seq_along(cl),function(i){
    co=cl[[i]]
    if (is.null(ini)) ini1=NULL else {
      x=ini[co]
      ini1=factor(x)
      ini1=as.double(ini1)
      names(ini1)=names(x)
    }
    if (length(co)>=s & max(l)>1) {
      y=treeCr(co,g,ini=ini1)
      return(y)
    }
    else return(co)
  })
  return(x)
}