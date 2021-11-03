adjL2Glite=function(AL){
require(igraph)
require(reshape2)
  
EL=melt(AL)
EL=as.matrix(EL)
G=graph_from_edgelist(EL, directed=F)
G=simplify(G)
#plot(G, vertex.size=1, main=paste("Graph at Cut Off",coff, collapse = " "))
return(G)
}