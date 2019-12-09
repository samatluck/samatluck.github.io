/*
 * ScalapackWrapper.h
 *
 *  Created on: Nov 13, 2013
 *      Author: fuji
 */

#ifndef SCALAPACKWRAPPER_H_
#define SCALAPACKWRAPPER_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#ifdef MKL_ILP64
#include <mkl_blas.h>
#include <mkl_pblas.h>
#include <mkl_blacs.h>
#include <mkl_lapack.h>
#include <mkl_scalapack.h>
#endif

class ScalapackWrapper {
public:
	ScalapackWrapper(int sizeOfMatrix,int sizeOfBlock, int sizeOfRowProcesses,int sizeOfColProcesses);
	virtual ~ScalapackWrapper();

	double getValue(int gi,int gj);
	bool setValue(int gi,int gj,double elm);

	double getValueRHS(int gi);
	bool setValueRHS(int gi,double elm);

	void generateTestSystem();
	void exportMat();

	bool isMyOwn(int gi,int gj);

	void pdgesv();
private:
	int sizeOfMatrix;
	int sizeOfBlock;
	int sizeOfRowProcesses;
	int sizeOfColProcesses;

	double *matA;
	double *vecB;
#ifdef MKL_ILP64
	MKL_INT myid,numproc;
	MKL_INT myrow, mycol;	
	MKL_INT *ippiv;
	MKL_INT descA[9], descB[9];
	MKL_INT nqrhs,np,nq;
	MKL_INT ictxt;
#else
	int myid,numproc;
	int myrow, mycol;
	int *ippiv;
	int descA[9], descB[9];
	int nqrhs,np,nq;
	int ictxt;
#endif
};

#endif /* SCALAPACKWRAPPER_H_ */
