/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use Lapack
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "stokeslet2dLapack.h"
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#ifdef MKL
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#else
extern "C" {
#include <atlas/cblas.h>
#include <atlas/clapack.h>
}
#endif
#endif

//Constructor
Rstokeslet2DLapack::Rstokeslet2DLapack(int numberOfParticles) :
				Rstokeslet2D(numberOfParticles) {
}

//Destructor
Rstokeslet2DLapack::~Rstokeslet2DLapack(){
}

/* Sovle Liear ststem */
/* Use Lapack */
bool Rstokeslet2DLapack::solve() {
	const int np = numberOfParticles;
	int info;

	/*solve Ax=b */
	/* A (amat) is n x n matrix */
	int nrhs = 1; /* # of RHS is 1 */
	int n = DIM * np;
	int lda = DIM * np;
	int ldb = DIM * np;
	/* The pivot indices that define the permutation matrix P;
	 row i of the matrix was interchanged with row IPIV(i). */
	int *ipiv = new int[DIM * np];
	
	/* set RHS */
	for (int i = 0; i < n; i++){
		foc[i] = vel[i];
	}
	
	/* solve with LAPACK */
#ifdef __APPLE__
	dgesv_(&n, &nrhs, mat, &lda, ipiv, foc, &ldb, &info);
#else
#ifdef MKL
	dgesv_(&n, &nrhs, mat, &lda, ipiv, foc, &ldb, &info);
#else
	info = clapack_dgesv(CblasColMajor, n, nrhs, mat, lda, ipiv, foc, ldb);
#endif
#endif

	/* check */
	if (info < 0) {
		std::cout << "the " << -info << "-th argument had an illegal value\n";
		return false;
	}
	if (info > 0) {
		std::cout << "U("<<info<<","<<info<<") is exactly zero.\n";
		return false;
	}

	/* free */
	delete[] ipiv;

	return true;
}

