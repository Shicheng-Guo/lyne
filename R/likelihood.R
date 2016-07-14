

mcTreeLike <- function(ali,edgeMat,transProb,eqFreq,takeLog=FALSE){
#=========================================================================

	nIntNodes = length(unique(edgeMat[,2]))
	nLeafs    = length(unique(edgeMat[!(edgeMat[,1] %in% edgeMat[,2]),1]))

	#- SOME CHECKS
	#=============

	#- Make sure edgeMat makes sense	
	if(sum(edgeMat[1:nLeafs,2] %in%  edgeMat[1:nLeafs,1]) > 0){
		stop("Error: illegal edgeMat\n")}
	for( i in (nLeafs+1):nrow(edgeMat)){
		if( ! edgeMat[i,1] %in% edgeMat[1:(i-1),2]  ) stop("Error: illegal edgeMat\n")
	}
	
	if(nrow(transProb) != ncol(transProb)) stop("transProb slices not square\n")
	if(nrow(transProb) != length(eqFreq))        stop("transProb to eqFreq mismatch\n")
	if(nrow(edgeMat) != dim(transProb)[3]) stop("edgeMat to transProb mismatch\n")
	if(length(unique(as.vector(edgeMat))) != nrow(ali)) stop("edgeMat to ali mismatch\n")

	#- GET THE LIKELIHOOD
	#====================
	
	ali[is.na(ali)] = -1 #- NAs to -1

	res = rcpp_likelihood(edgeMat,nIntNodes,ali,dim(transProb)[1],transProb,eqFreq,takeLog)
	return(as.vector(res)) 

}

mcTreeLike.3states.GTR.bhom <- function(pars, ali,mult=1,edgeMat,takeLog=FALSE,logPars=FALSE,sumUp=TRUE){
#========================================================================================================
	
	#- FOR BRANCH-HOMOGENUOUS (same Q) GTR using 3 STATES
	#====================================================

	#pars[1:3]            #- equilibrium frequencies
	#pars[4:6]            #- the alphas: 12,13 and 23
	#pars[7:length(pars)] #- the branch lengths

	if(logPars==TRUE) pars = exp(pars)


	#- 0. Normalize frequencies and alphas (easier for optim)
	pars[1:3] = pars[1:3]/sum(pars[1:3])
	pars[4:6] = pars[4:6]/sum(pars[4:6])
	
	#- 1. Make the Q
	QEig      = pars2QEig.3states.GTR(pars[1:6])

	#- 2. Make the transition probabilities
	transProb = aperm(maply(pars[-(1:6)],function(x) QEig2P(QEig,x)),perm=c(2,3,1))
	
	#- 3. Calculate the likelihood for each column
	likes     = mcTreeLike(ali,edgeMat,transProb,pars[1:3],takeLog=takeLog)
	
	#- 3.1 multiplicities
	likes     = likes*mult

	#- 4. Sum up and be done
	res = sum(likes)
	if(!sumUp) return(likes) #- returns log-likelihoods per observation
	return(res)              #- default returns summed-up log likelihood
}

mcTreeLike.3states.inv <- function(ali,mult=1,sumUp=TRUE,multMult=TRUE){
#========================================================================

	#- Just use empirical frequencies

	if(length(mult)==1) mult = rep(mult,ncol(ali)) #- FIXME

	all.1.ind = which(apply(ali,2,function(x) all(x==1)))
        all.2.ind = which(apply(ali,2,function(x) all(x==2)))
        all.3.ind = which(apply(ali,2,function(x) all(x==3)))

	n.1   = 1 ; if(length(mult)>1) n.1 = sum(mult[all.1.ind])
	n.2   = 1 ; if(length(mult)>1) n.2 = sum(mult[all.2.ind])
	n.3   = 1 ; if(length(mult)>1) n.3 = sum(mult[all.3.ind])
	n.inv = length(all.1.ind)*n.1 + 
		length(all.2.ind)*n.2 + 
		length(all.3.ind)*n.3
	p.1   = length(all.1.ind)/n.inv*n.1
	p.2   = length(all.2.ind)/n.inv*n.2
	p.3   = length(all.3.ind)/n.inv*n.3

	ind.1 = which(apply(ali,2,function(x) all(x[!is.na(x)]==1)))
	ind.2 = which(apply(ali,2,function(x) all(x[!is.na(x)]==2)))
	ind.3 = which(apply(ali,2,function(x) all(x[!is.na(x)]==3)))
	
	colLik        = rep(-Inf,ncol(ali))
	colLik[ind.1] = log(p.1)
	colLik[ind.2] = log(p.2)
	colLik[ind.3] = log(p.3)
	if(multMult) colLik        = colLik * mult

	if(!sumUp) return(colLik)
	return(sum(colLik)) #- will be -Inf, except for only inv sites
}

gamma.make.rates <- function(alpha,ncat=3){
#==========================================

  	shape     = alpha   #- alpha is the shape parameter
        scale     = 1/alpha #- so we have mean one for the distribution
        midcuts   = qgamma(0:ncat/ncat,shape=shape,scale=scale)[-c(1,ncat+1)]

        #- need the expectation in each of the ncat bins
        starts = c(0,midcuts)
        ends   = c(midcuts,Inf)
        rates  = rep(NA,ncat)
        ifu    = function(x) x*dgamma(x,shape=shape,scale=scale)
        for(i in 1:ncat){
                tmp = integrate(ifu,starts[i],ends[i])
                if(tmp$message=="OK") rates[i] = tmp$value
        }
        rates = rates * ncat/sum(rates) #- so that mean(rates) = 1 for sure
        return(rates)

}

mcTreeLike.3states.GTR.bhom.gamma <- function(pars,ali,mult=1,edgeMat,takeLog=FALSE,logPars=FALSE,sumUp=TRUE,details=FALSE){
#===========================================================================================================================

	if(logPars == TRUE) pars = exp(pars) ; pars[1] = round(pars[1]) #- number of cats needs to be integer

	rates        = gamma.make.rates(pars[2],pars[1])
	#- BELOW: mult *must not* be passed here, needs to be set to 1, this is not a default argument
	colLikByRate = maply(rates,function(x) mcTreeLike.3states.GTR.bhom(pars[-c(1,2)]*c(rep(1,6),rep(x,length(pars)-8)),ali,mult=1,edgeMat,sumUp=FALSE,logPars=FALSE))	
	colLik       = log(colMeans(exp(colLikByRate)))	
	colLik	     = colLik*mult
	if(!sumUp) {
		if(!details) return(colLik)
		return(colLikByRate)
	}
	return(sum(colLik))

}

mcTreeLike.3states.GTR.bhom.gamma.inv <- function(pars,ali,mult=1,edgeMat,takeLog=FALSE,logPars=FALSE,sumUp=TRUE,details=FALSE){
#===============================================================================================================================

	if(logPars == TRUE) pars = exp(pars) ; pars[2] = round(pars[2]) #- number of rate cats needs to be integer

	if(pars[1] > 1) stop("illegal mixture coefficient: too much inv sites.")
	if(pars[1] < 0) stop("illegal mixture coefficient: too little inv sites.")

	#- the likelihood for non-invariable sites 
	rates        = gamma.make.rates(alpha=pars[3],ncat=pars[2])
	#- BELOW: mult *must not* be passed here, needs to be set to 1, this is not a default argument
	colLikByRate = maply(rates,function(x) mcTreeLike.3states.GTR.bhom(	pars[-c(1,2,3)]*c(rep(1,6),rep(x,length(pars)-9)),
										ali,mult=1,edgeMat,sumUp=FALSE,logPars=FALSE))	
	colLik       = colMeans(exp(colLikByRate))
	
	#- the likelihood for the invariable sites 
	colLik.inv   = exp(mcTreeLike.3states.inv(ali,mult=mult,sumUp=FALSE,multMult=FALSE))

	#- the mixture likelihood
	colLik.res   = pars[1]*colLik.inv + (1-pars[1])*colLik
	#colLik.res[colLik.inv == 0] = colLik[colLik.inv==0]
	
	#- be done
	colLik.res   = log(colLik.res)*mult

	if(!sumUp){
		if(!details) return((colLik.res))
		res = (rbind(log(colLik.inv),colLikByRate))
		rownames(res) = c(0,1:nrow(colLikByRate))
		return(res)
	}
	return(sum((colLik.res)))

}




optimize.mcTreeLike.3states.GTR.bhom <- function(pars,ali,mult=1,edgeMat,takeLog=FALSE){
#=======================================================================================

	res       = optim(	log(pars)					,
				fn 	= mcTreeLike.3states.GTR.bhom		,
				method 	="L-BFGS-B"				,
				lower 	= rep(-12,length(pars))			,
				upper  	= rep(6,length(pars))			,
				control = list(fnscale=-1,factr=1e12,pgtol=1E-3),
				ali     = ali					,
				mult    = mult					,
				edgeMat = edgeMat				,
				takeLog = takeLog				,
				logPars = TRUE)
	return(res)
}

optimize.mcTreeLike.3states.GTR.bhom.gamma <- function(pars,ali,mult=1,edgeMat,takeLog=FALSE){
#=============================================================================================

	ncats = pars[1]
	ofu   <- function(pars,ali,mult,edgeMat) {
		#- not optimize over ncats ; look up from parent env.
		res = mcTreeLike.3states.GTR.bhom.gamma(pars=c(log(ncats),pars),ali=ali,mult=mult,edgeMat=edgeMat,logPars=TRUE,takeLog=TRUE)
		sum(res) }
  	res       = optim(      log(pars[-1])                                 	,
                                fn       = ofu           			,
                                method  ="L-BFGS-B"                             ,
                                lower   = c(-3,rep(-12,length(pars-1)))        	,
                                upper   = rep(6,length(pars))           	,
                                control = list(fnscale=-1,factr=1e12,pgtol=1E-3),
                                ali     = ali                                   ,
                                mult    = mult                                  ,
                                edgeMat = edgeMat                               )
        return(res)

}


optimize.mcTreeLike.3states.GTR.bhom.gamma.inv <- function(pars,ali,mult=1,edgeMat,takeLog=FALSE){
#=================================================================================================

	ncats = pars[2]
	ofu   <- function(pars,ali,mult,edgeMat) {
		print(pars)
		res = mcTreeLike.3states.GTR.bhom.gamma.inv(pars=c(pars[1],log(ncats),pars[-1]),ali=ali,mult=mult,edgeMat=edgeMat,logPars=TRUE,takeLog=FALSE)
		print(res)
		sum(res) }
  	res       = optim(      log(pars[-2])                                   ,
                                fn      = ofu           			,
                                method  ="L-BFGS-B"                             ,
                                lower   = c(-3,-3,rep(-12,length(pars-3)))      ,
                                upper   = c(-1E-10,rep(6,length(pars-2)))       ,
                                control = list(fnscale=-1,factr=1e12,pgtol=1E-3),
                                ali     = ali                                   ,
                                mult    = mult                                  ,
                                edgeMat = edgeMat                               )
        return(res)

}




