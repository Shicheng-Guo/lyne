// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace Rcpp; 

arma::vec logProd(arma::vec,arma::mat);

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// ######################
// LIKILIHOOD CALCULATION
// ###################### 
// ################################################################
// [[Rcpp::export]]
arma::vec rcpp_likelihood(	const	arma::imat mEM		,
				const	int nIntNodes		,
				const	arma::imat mAli		,
				const	int nStates		,
					NumericVector vProb	,
				const	arma::rowvec vPi	,
				const 	bool takeLog		){
// ###############################################################

// * mEM is edge matrix in bottom-up order leaf->root
// * mAli are lignment columns with states {-1,0,1,2,...,nStates-1}
//   	mAli has -1 for missing data
//	mAlis and mEM need to be aligned: entries in mEM refer to rows in mAli
//	
// * aliCol has -1 for missing/unobserved data ; FIXME
// * vProb matches the order in mEM:
// 	vProb is from a 3-dimensional R-arrray: nStates x nStates x nNodes
// 	the nNodes slices get indexed by mEM
// * leaf-vertex index <= nLeafs; i.e.:
//	max(mEM[1:nLeafs,1] <= nLeafs) the from-indices are small, only from-leafs
// 	min(mEM[1:nLeafs,2]  > nLeafs) the to-indices are large, no to-leafs
// * takeLog == 1 takes log of partial likelihoods on tree. SLOWER, no mat-mult.

	int nNodes   = mAli.n_rows; 		// number of nodes in the tree. All could be observed. 
	int nAliCols = mAli.n_cols; 		// number of alignment columns
	int nEdges   = mEM.n_rows ; 		// number of edges in tree
	int nLeafs   = nNodes - nIntNodes ; 	// root is "internal" node here
	
	// re-shape the transition matrices
	//---------------------------------
	arma::cube cProb(vProb.begin(), nStates, nStates, nEdges, false); // 3D array
	arma::cube lcProb;
	if(takeLog==true){ // make log transition probs
		lcProb = arma::trunc_log(cProb); 
	}
	
	// initialize partial likelihoods
	//-------------------------------
	arma::mat pLike = arma::ones<arma::mat>(nNodes, nStates);
	arma::vec rLike = arma::zeros<arma::vec>(nAliCols);
	arma::mat pLiket;

	//- for rach alignment column
	//---------------------------
	for(int j=0; j<nAliCols; j++){
	
		if(j>0){
			pLike.ones();
		}
	
		//- partial likelihood for all nodes 
		//---------------------------------
		for (int i=0; i < nNodes; i++){
			pLike.row(i).zeros() ; 			// nothing goes by default  
			if (mAli(i,j) == -1){  
	      			pLike.row(i).ones(); 		// -1 denotes a missing value, anything goes
			} else {
				pLike(i, mAli(i,j)-1) = 1.0; 	// clamp to observed value, only those go, -1 because zero-based here
			}
		}
	
		if(takeLog == true ){
			pLike  = arma::trunc_log(pLike); //- FIXME: bad zero...
			pLiket = pLike.t();        //- save some transposing downstream 
		}
	
		// partial likelihoods for internal nodes
		//---------------------------------------
		for (int i=0; i < nEdges; i++){
			if(takeLog == false ){
				// BELOW: FIXME TRANSPOSE
				pLike.row( mEM(i,1)-1) %= pLike.row(mEM(i,0)-1) * cProb.slice(i).t(); // * is matrix-mult, % element-wise	
				// the -1 in  mEM(i,1)-1 is because mEM is assumed to contain indices 1...nNodes, but here we are zero-based
			} else {
				pLiket.col( mEM(i,1)-1) += logProd(pLiket.col(mEM(i,0)-1),lcProb.slice(i).t());
				// FIXME: SLOW!
			}
		}
	
		// save likelihood for this alignment column
		if(takeLog == true ){
			rLike(j) = arma::trunc_log(sum(arma::trunc_exp(arma::trunc_log(vPi.t()) + pLiket.col(mEM(nEdges-1,1)-1)))); // mEM(nEdges-1,1) "points" to root  
		} else {
			rLike(j) = sum(vPi % pLike.row(mEM(nEdges-1,1)-1)); // mEM(nEdges-1,1) "points" to root
			rLike(j) = arma::trunc_log(rLike(j));
		}

	}
	return(rLike);
}
	
//########################################
arma::vec logProd(	arma::vec lpLik, 
			arma::mat  lP  ){
//#######################################	

		arma::vec res(lpLik.n_elem);
			
		for(int i=0; i<lP.n_cols; i++){
			res(i) = sum(arma::trunc_exp(lpLik  + lP.col(i)));
		}
		
		return(arma::trunc_log(res));
}

//rEM  = rbind( c(1,5),c(2,5),c(3,6),c(4,6),c(5,7),c(6,7))
//rAli = matrix(sample(0:3,1000000,rep=T),nrow=4)
//rAli = rbind(rAli,matrix(-1,ncol=ncol(rAli),nrow=3))
//pm   = matrix(1/4,ncol=4,nrow=4)
//PM   = array(pm,dim=c(4,4,nrow(rEM)))
//rPi  = rep(1/4,4)
//tst  = rcpp_likelihood(rEM,3,rAli,4,PM,rPi,0)



