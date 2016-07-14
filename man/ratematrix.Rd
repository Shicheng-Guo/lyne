
\name{ratematrix}
\alias{pars2Q.3states.GTR}
\alias{pars2QEig.3states.GTR}
\alias{Q2QEig}
\alias{QEig2Q}
\alias{QEig2P}
\alias{Q2P}

\title{
	Translate between parameters, rate matrices and transition probilities. 
}
\description{
Calculate ratematrix from GTR parameterization for 3 states, get Eigen decomposition and
calculate transition probabilities.
}

\usage{
pars2Q.3states.GTR(pars)
pars2QEig.3states.GTR(pars)
Q2QEig(Q)
QEig2Q(QEig)
QEig2P(QEig,t)
Q2P(Q,t)
}

\arguments{
  \item{pars}{
	Vector: Parameters describing the substitution model for each edge on the tree.
  }
  \item{Q}{
	Matrix: Ratematrix.
  }
  \item{QEig}{
	List: Three components: \code{d}: Eigenvalues, \code{Ui} Matrix of Eigenvectors, 
	\code{U} Inverse of \code{Ui}
  }
  \item{t}{
	Numeric: Time, i.e. the expected number of substitutions along that branch in the tree.
 }
}
\details{
Parameters \code{pars} are expanded to a rate matrix by using the first three values as equilibrium
frequencies, and the following three values as rate parameters (1->2, 1->3 and 2->3).
}

\value{
\code{pars2Q.3states.GTR} Returns a 3x3 rate matrix.
\code{pars2QEig.3states.GTR} Returns the Eigen decomposition of a  3x3 rate matrix.
\code{Q2QEig} Returns the Eigen decomposition of Q.
\code{QEig2Q} Returns the rate matrix corresponding to the Eigen decomposition \code{QEig}.
\code{QEig2P} Returns the transition matrix.
\code{Q2P} Returns the transition matrix.

 
}
