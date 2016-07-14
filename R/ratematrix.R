
pars2Q.3states.GTR <- function(pars){
#===================================

#- Make a (NOT eigen-decomoposed) rate matrix from six parameters:	
	
	#- three eq-freq pi1,pi2,pi3
        #- three alphas  al12,al13,al23
	
	Qm = rbind( c(0,              pars[2]*pars[4],pars[3]*pars[5]),
                    c(pars[1]*pars[4],0,              pars[3]*pars[6]),
                    c(pars[1]*pars[5],pars[2]*pars[6],0              ))

        diag(Qm) = -rowSums(Qm)
        Qm       = Qm/( (-1)*sum(diag(Qm)*pars[1:3]) ) #- eq subs rate one

	return(Qm)
}

pars2QEig.3states.GTR <- function(pars){
#======================================

#- Make a (eigen-decomoposed) rate matrix from six parameters:
	
	#- three eq-freq pi1,pi2,pi3 
	#- three alphas  al12,al13,al23

	Qm  = pars2Q.3states.GTR(pars)
	tmp = eigen(Qm)

	Q = list()
	Q$d  = tmp$values 
	Q$Ui = tmp$vectors
	Q$U  = solve(tmp$vectors)

	return(Q)
}

Q2QEig <- function(Q){
#=====================

#- Eigen-decompose rate matrix

	tmp     = eigen(Q)
	QEig    = list()
        QEig$d  = tmp$values
        QEig$Ui = tmp$vectors
        QEig$U  = solve(tmp$vectors)

        return(QEig)
}

QEig2Q <- function(QEig){
#========================

#- Get Q from its eigen decomposition

	Q = QEig$Ui %*%diag(QEig$d) %*% QEig$U
	return(Q)
}


QEig2P <- function(QEig,t){
#==========================

#- Make a transition matrix from an (eigen-decomposed) rate matrix

	P = zapsmall( QEig$Ui %*% diag(exp(QEig$d*t)) %*% QEig$U )
	P = P/rowSums(P) #- FIXME: fix rounding error
	return(P)
}

Q2P <- function(Q,t){
#=====================

#- Make a transition matrix from an (NOT eigen-decomposed) rate matrix
	
	QEig = Q2QEig(Q) 
        P    = zapsmall( QEig$Ui %*% diag(exp(QEig$d*t)) %*% QEig$U )
	P = P/rowSums(P) #- FIXME: fix rounding error
	return(P)
}


