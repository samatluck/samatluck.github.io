/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1

 Use Eigen3
*/
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>

#define DIM 2
#define EPSILON 0.005 /* blob size*/

class Rstokeslet2DEigen3{
public:
	Rstokeslet2DEigen3(int numberOfParticles);
	virtual ~Rstokeslet2DEigen3();
	
	void mkMatrix();
	bool solve();
	void getVelocities(int, double *iLoc, double *iVel);

	double *getLocationArray();
	double *getVelocityArray();
	double *getForceArray();
	double *getMatrix();

private:	
	double term1(double, double);
	double term2(double, double);

	int numberOfParticles; // number of particlesed:

	void matZero(); // zeros matrix
	void setElmMat(int i,int j,double elm); // insert a value at (i,j)
	void addElmMat(int i,int j,double elm); // add a value at (i,j)
	
	Eigen::VectorXd loc;
	Eigen::VectorXd vel;
	Eigen::VectorXd foc;
	Eigen::MatrixXd mat;
};
