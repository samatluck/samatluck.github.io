/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1

 Use Scalapack
*/
#include <iostream>
#include "ScalapackWrapper.h"

#define DIM 2
#define EPSILON 0.005 /* blob size*/
#define BLOCK_SIZE 64

class Rstokeslet2DScalapack{
public:
	Rstokeslet2DScalapack(int numberOfParticles,int nump,int numq);
	virtual ~Rstokeslet2DScalapack();
	
	void mkMatrix();
	bool solve();
	void getVelocities(int, double *iLoc, double *iVel);

	double *getLocationArray();
	double *getVelocityArray();
	double *getForceArray();

	void collectResults();
private:	
	double term1(double, double);
	double term2(double, double);

	int numberOfParticles; // number of particlesed:

	void matZero(); // zeros matrix
	void setElmMat(int i,int j,double elm); // insert a value at (i,j)
	void addElmMat(int i,int j,double elm); // add a value at (i,j)
	
	double *loc;/* particle location */
	double *vel;/* particle velocity */
	double *foc;/* particle force */

	ScalapackWrapper sclwrap;
	
	int myid;
	int numproc;
};
