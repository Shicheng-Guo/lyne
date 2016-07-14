getEdgeMatrix<- function(graph,edges=FALSE){
#===========================================

#- Make a tree traversal to calc the likelihood a la felsenstein


        #- root and leafs of the graph
        root     = V(graph)[which(degree(graph,mode="in") == 0)]
        leafs    = V(graph)[which(degree(graph,mode="out") == 0)]
        leaf.ind = which(V(graph) %in% leafs)

        #- get valid vertex order
        vert.ord     = graph.dfs(graph,root)$order
        vert.ord     = rev(vert.ord) #- reverse dfs
        vert.ind.ord = vert.ord[!vert.ord %in% leaf.ind]

        #- get the edges for each vertex we'll visit
        E.ord.list   = lapply(vert.ind.ord, function(x) E(graph)[from(x)])
        edge.ind.ord = unlist(E.ord.list) #- indexes of the edges in order to visit

        #- for each edge, get the vertices in the right order
        lefu <- function(x){
                from.ind  = as.numeric(V(graph)[from(x)]) #- one number
                to.ind    = as.numeric(V(graph)[to(x)]) #- one (or many) numbers
                return(cbind(from.ind,to.ind))
        }
        D.ord.list = lapply(E.ord.list,lefu)

        #- put everything in a matrix and be done
        Dmat = matrix(unlist(lapply(D.ord.list,function(x) t(x))),ncol=2,byrow=T)
        rMat =  cbind(edge.ind.ord,Dmat)
	
	#- reverse edges so they point to root
	rMat = rMat[,c(1,3,2)]

	colnames(rMat) = c("edge","node.from","node.to")

	if(!edges) rMat = rMat[,-1]

        return(list(traversal = rMat, leafs=leaf.ind))
}


