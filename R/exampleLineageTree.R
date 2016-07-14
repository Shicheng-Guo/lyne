
exampleLineageTree<- function(){
#===============================

#- This function makes an  igraph object representing the Bock et al. 2012 blood cell lineages
#  (http://dx.doi.org/10.1016/j.molcel.2012.06.019)

	el = rbind(	c("hsc","mpp1"),
			c("mpp1","mpp2"),
			c("mpp2","cmp"),
			c("cmp","gmp"),
			c("gmp","mono"),
			c("gmp","gran"),
			c("cmp","mep"),
			c("mep","ery"),
			c("mpp2","clp"),
			c("clp","cd4"),
			c("clp","cd8"),
			c("clp","bcell"))

	g = graph.edgelist(el)

	layout = rbind(	c(1.25,6/6)	,
			c(1.25,5/6)	,
			c(1.25,4/6)	,
			c(1.875,3/6)	,
			c(2.25,2/6)	,
			c(2.50,1/6)	,
			c(2.00,1/6)	,
			c(1.50,2/6)	,
			c(1.50,1/6)	,
			c(0.50,2.5/6)	,
			c(0.00,1/6)	,
			c(0.50,1/6)	,
			c(1.00,1/6)	
			)

	#- graph is the graph, layout is for plotting
	return(list(graph=g,layout=layout))
}

