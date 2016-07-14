
\name{mcTreeLike}
\alias{mcTreeLike}
\alias{mcTreeLike.3states.inv}
\alias{mcTreeLike.3states.GTR.bhom}
\alias{mcTreeLike.3states.GTR.bhom.gamma}
\alias{mcTreeLike.3states.GTR.bhom.gamma.inv}
\alias{mcTreeLike.3states.GTR.bhom.gamma.inv}
\alias{optimize.mcTreeLike.3states.GTR.bhom}
\alias{optimize.mcTreeLike.3states.GTR.bhom.gamma}
\alias{optimize.mcTreeLike.3states.GTR.bhom.gamma.inv}

\title{
	Calculate the likelihood for Markov chain models on trees. 
}
\description{
Likelihood for a given alignment and tree model (tree and transition probabilities between states for 
each edge). Sums over missing data (elimination algorithm). General trees, any node 
can be missing.
}

\usage{
mcTreeLike(ali, edgeMat, transProb, eqFreq, takeLog=FALSE)
mcTreeLike.3states.inv(ali, mult=1, sumUp=TRUE, multMult=TRUE)
mcTreeLike.3states.GTR.bhom(pars, ali, mult=1, edgeMat, takeLog=FALSE, logPars=FALSE, sumUp=TRUE)
mcTreeLike.3states.GTR.bhom.gamma(pars, ali, mult=1, edgeMat, takeLog=FALSE, logPars=FALSE, sumUp=TRUE, details=FALSE)
mcTreeLike.3states.GTR.bhom.gamma.inv(pars, ali, mult=1, edgeMat, takeLog=FALSE, logPars=FALSE, sumUp=TRUE, details=FALSE)
optimize.mcTreeLike.3states.GTR.bhom(pars,ali,mult=1,edgeMat,takeLog=FALSE)
optimize.mcTreeLike.3states.GTR.bhom.gamma(pars,ali,mult=1,edgeMat,takeLog=FALSE)
optimize.mcTreeLike.3states.GTR.bhom.gamma.inv(pars,ali,mult=1,edgeMat,takeLog=FALSE)

}

\arguments{
  \item{pars}{
	Vector: Parameters describing the substitution model on the tree (see Details)
  }
  \item{ali}{
	Matrix: (number of tree nodes) x (number of observations). Contains integers, each representing a state in the Markov chain.
	Columns are observations, rows tree-nodes. 	
  }
  \item{mult}{
	Vector: Multiplicity of alignment columns. Allows the use of sufficient statistics.
  }
  \item{edgeMat}{
	Matrix: integer edge matrix of the tree in bottom-up order. Two columns: first is the from-nodes (indexes), second the to-nodes for
	each edge. Tree needs to be rooted. Nodes are integers that correspond to rows \code{ali}.
  }
  \item{transProb}{
	Array: transition matrix for each edge. (number of states) x (number of states) x (number of edges in tree)
  }
  \item{eqFreq}{
	Vector: The equilibrium frequencies at the root node of the tree. 
  }
  \item{takeLog}{
	Logical: Should log be taken while calculating likelihood?
  }
  \item{logPars}{
	Logical: Are parameters \code{pars} provided on \code{log}-scale?
  }
  \item{sumUp}{
	Logical: Should the likelihood be summed over alignment columns, or should the likelihood for each column be reported separately?
  } 
  \item{details}{
	Logical: For mixture models: Should the likelihood of each component be returned, or should the components be (weighted and) summed?
  }
  \item{multMult}{
	Logical: If multiplicities are given (\code{mult != 1}), and if \code{sumUp != TRUE}, should the log-likelihood for each observation be 
	multiplied by its multiplicity?  
  }	
}

\details{
	All function calculate the log-likelihood, given a model specification. \code{mcTreeLike} is the most generic, and is used by all the 
	other functions (except \code{mcTreeLike.3state.inv}) after converting parameter values into transition probabilities. The matrix 
	\code{edgeMat} contains integers referencing rows in \code{ali} as the tree nodes. It needs to be provided in bottom-up order, such that
	its rows provide a traversal from the leafs to the root of the tree.\cr

	\code{cTreeLike.3state.inv} Calculates the likelihood assuming invariable states only. It uses the frequencies of completely-observed
	(i.e., no missing data )invariable states in \code{ali} and does not need additional parameters. 
	Likelihood of observations with more than one state is zero. 
	Observations with only one state and missing data are assigned the likelihood of the compatible invariant state.\cr

	\code{mcTreeLike.3state.GTR.bhom} Branch-homogeneous model. The parameter vector is organized as follows:\cr
	\code{pars[1:3]}: Equilibrium frequencies\cr
	\code{pars[4:6]}: Rate parameters (1->3, 1->3, 2->3) \cr
	\code{pars[-(1:6)]}: Branch lengths of the tree\cr
	\cr
	
	The functions with the \code{optimize} prefeix optimize the respective likelihood function, given data and initial parameters \code{pars}.
}

\value{
\code{treeLike} returns a vector containing the likelihood for each column in the alignment matrix \code{ali}. If \code{takeLog == TRUE} it is on log scale.	
}
