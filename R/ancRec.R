
recCol <- function(col,state,pars,rates=NULL,aligned=FALSE){
#===========================================================
	
	#- PARS not on log
	#- gamma.inv pars expected

	if(all(is.na(col)))  return(list(rec=col,lik=NA))
        if(!any(is.na(col))) return(list(rec=col,lik=NA))

        #- invariant site
        if(state==0){
                col[is.na(col)] = unique(col[!is.na(col)])
                return(list(rec=col,lik=NA))
        }
	
	#- variable site
	LT        =  exampleLineageTree()
	tree      = LT$graph
	traversal = getEdgeMatrix(tree)
	pars[-(1:8)] = pars[-(1:8)][order(traversal$traversal[,1])] #- OLD PARS NEEDED
	traversal$traversal = traversal$traversal[,c(1,3,2)]        #- OLD edgeMat needed
	colnames(traversal$traversal)  = colnames(traversal$traversal)[c(1,3,2)]


	if(is.null(rates)){
		rates = gamma.make.rates(pars[2],3)
	}
	pars  = pars[-(1:2)]
	pars[-(1:6)] = pars[-(1:6)] * rates[state]

	bg     = pars[1:3]/sum(pars[1:3])
        al     = pars[4:6]/sum(pars[4:6])
        ts     = pars[-(1:6)]
	QEig   = pars2QEig.3states.GTR(c(bg,al))
	Ps     = lapply(ts,function(x) QEig2P(QEig,x))
	logPts = lapply(Ps,function(x)log(t(x)))

	rec = logLikeAndRec(    col,
                                                        traversal$traversal,
                                                        leafs=traversal$leafs,
                                                        liks=NULL,
                                                        recs=NULL,
                                                        logPts,
                                                        bg)

        return(rec)
}

logLikeAndRec <- function(dat,traversal,leafs=NULL,liks=NULL,recs=NULL,logPts,root.eq){
#=======================================================================================

        #- dat:         The data. Needs to be aligned such that dat[i] ~ V(mod)[i].
	#		Missing values will be reconstructed (joint maximum likelihood)
        #- traversal:   Tree traversal from leafs to root. Integer matrix with 3 columns:
        #               edge: The edge (needs to match E(mod)) for which we are doing the likelihood
        #               node.from: the node from which the edge is coming (index thereof)
        #               node.to: the node to which the edge is going (index thereof)
        #- leafs:       Index of the leafs in the tree ( leafs[i] ~ V(mod)[i] )
 	#- liks:        Initial values for the likelihoods on the dag. Should be NULL.
	#- recs:	Initial values for the reconstructions on the dat. Shold be NULL>
        #- logPs:       The *transposed* log transition matrixes; logPts[[i]] ~ E(mod)[i]
        #- root.eq:     The equilibrium frequencies at the root

        #- do we need to calculate the leafs?
        if(is.null(leafs)) leafs = traversal[,"node.to"][!traversal[,"node.to"] %in% traversal[,"node.from"]]

	#- initialize the reconstructions 
	leafs.obs = leafs[!is.na(dat[leafs])]
	leafs.mss = leafs[ is.na(dat[leafs])]

	dat.leafs.obs = dat[leafs.obs]
	dat.leafs.mss = dat[leafs.mss]
	
	#- initialize. We'll use PUPKO et AL., MBE 2000 (with according generalizations) for joint ML reconstruction
	if(is.null(recs)){
                recs  = matrix(NA,ncol=3,nrow=nrow(traversal)+1)
       	 }
	
	if(is.null(liks)){
		liks  = matrix(-Inf,ncol=3,nrow=nrow(traversal)+1)
        }

	#- fill up non-root likelihoods and reconstructions
	for(i in 1:nrow(traversal)){
		n.to        = traversal[i,"node.to"]
                n.fo        = traversal[i,"node.from"]
		P.fo        = logPts[[traversal[i,"edge"]]]
	
		if(n.to %in% leafs.obs){
			if(is.na(dat[n.fo])){
				recs[n.to,] = rep(dat[n.to],ncol(recs)) #- always reconsturct to observation
				liks[n.to,] = P.fo[,dat[n.to]]          #- from anything to the observation
			} else {
				recs[n.to,dat[n.fo]] = dat[n.to]                 #- reconsturct to observation from observed parent
                                liks[n.to,]          = rep(-Inf,ncol(recs))      #- from parent obervation to child observation
				liks[n.to,dat[n.fo]] = P.fo[dat[n.fo],dat[n.to]]
			}
		} else if(n.to %in% leafs.mss){
			if(is.na(dat[n.fo])){
				recs[n.to,] = apply(P.fo,1,which.max)   #- the state at the leaf with the maximum like for each parent 
				liks[n.to,] = apply(P.fo,1,max)	        #- the according likelihood
			} else {
				recs[n.to,dat[n.fo]] = which.max(P.fo[dat[n.fo],])
				liks[n.to,]          = rep(-Inf,ncol(recs))
				liks[n.to,dat[n.fo]] = max(P.fo[dat[n.fo],])
			}
		} else {
			#= the children of the to node
			n.kids = traversal[,"node.to"][traversal[,"node.from"]==n.to]
			l.kids = liks[n.kids,,drop=FALSE]
						
			if(is.na(dat[n.to])){ #- current node not observed
				if(is.na(dat[n.fo])){
					P.tmp = P.fo
					for(j in 1:nrow(l.kids)) P.tmp = t(t(P.tmp) + l.kids[j,]) 

					recs[n.to,] = apply(P.tmp,1,which.max)
					liks[n.to,] = apply(P.tmp,1,max)
				} else {
					liks[n.to,dat[n.fo]] = max(P.fo[dat[n.fo],] + colSums(l.kids))
					recs[n.to,dat[n.fo]] = which.max(P.fo[dat[n.fo],] + colSums(l.kids))
				}			
			} else { #- current node oberved
				if(is.na(dat[n.fo])){ #- parent node not observed
					recs[n.to,] = rep(dat[n.to],ncol(recs)) #- reconstruct to observed value
					liks[n.to,] = P.fo[,dat[n.to]] + sum(l.kids[,dat[n.to]])
					
                                } else { #- parent node observed
					recs[n.to,dat[n.fo]] = dat[n.to]
					liks[n.to,]          = rep(-Inf,ncol(recs))
					liks[n.to,dat[n.fo]] = P.fo[dat[n.fo],dat[n.to]] + sum(l.kids[,dat[n.to]])
                                }
			}
		}
	}

	#- reconstruct the root (we are node.from)
	i      = nrow(traversal)
	n.to   = traversal[i,"node.from"]
	n.kids = traversal[,"node.to"][traversal[,"node.from"]==n.to]
        l.kids = liks[n.kids,,drop=FALSE]
	
	if(is.na(dat[n.to])){
		recs[n.to,] = rep(which.max(colSums(l.kids)+log(root.eq)) ,ncol(recs))
	} else {
		recs[n.to,] = rep(dat[n.fo],ncol(recs)) #- reconstruct root to observations 
	}
	
	#- the likelihood of the reconstruction
	rec.lik = (colSums(l.kids)+log(root.eq))[recs[n.to,1]]	

	#- traverse agin (other direction) to reconstruct
	for(i in nrow(traversal):1){
		n.to         = traversal[i,"node.to"]
                n.fo         = traversal[i,"node.from"]
		recs[n.to,]  = rep(recs[n.to,recs[n.fo,1]],3) #- recs[n.fo,] has already been reconstructed here
	}

	#- our reconstruction
	rec.dat = apply(recs,1,unique)

	return(list(rec=rec.dat,lik=rec.lik))
}




















