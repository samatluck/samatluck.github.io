/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1

 Use PetSc
*/
#include <iostream>
#include <petsc.h>

#define DIM 2
#define EPSILON 0.005 /* blob size*/

class Rstokeslet2DPetSc{
public:
	Rstokeslet2DEigen3(int numberOfParticles);
	virtual ~Rstokeslet2DEigen3();
		
	void mkMatrix();
	bool solve();
	void getVelocities(int, double *iLoc, double *iVel);

	double *getLocationArray();
	double *getVelocityArray();
	double *getForceArray();
	
	void restoreLocationArray();
	void restoreVelocityArray();
	void restoreForceArray();

	int myStart;
	int myEnd;
	int mySize;
private:	
	double term1(double, double);
	double term2(double, double);

	int numberOfParticles; // number of particlesed:

	void matZero(); // zeros matrix
	void setElmMat(int i,int j,double elm); // insert a value at (i,j)
	void addElmMat(int i,int j,double elm); // add a value at (i,j)
	
	Vec loc;
	Vec vel;
	Vecfoc;
	Mat mat;
	
	int myid;
	int numproc;
};
