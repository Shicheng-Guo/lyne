
\name{tree}
\alias{tree}
\alias{getEdgeMatrix}
\alias{exampleLineageTree}

\title{
	Tree-related functions
}
\description{
Get an example of a lineage tree. Given a lineage tree, get an edge matrix accordingly
in bottom-up order (from leaves to root).
}


\usage{
exampleLineageTree()
getEdgeMatrix(graph,edges=FALSE)
}

\arguments{
  \item{graph}{
  	igraph-class: Object containing the lineage tree. Needs to be a rooted tree.
  }
  \item{edges}{
	Logical: Should the result contain a column with the edge index in \code{tree}
	corresponding to the respective row?
  }
}
\details{
\code{exampleLineageTree} returns the hematopoietic cell lineage tree from Bock et al. (2012). 
}
\value{
\code{getEdgeMatrix} Edge matrix in bottom-up order.\cr
\code{exampleLineageTree} List with two components: \code{traversal}: igraph-class object containing cell lineage tree.
  \code{layout}: layout for plotting the tree.

 
}
