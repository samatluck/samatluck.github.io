/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 */
#define DIM 2
#define EPSILON 0.005 /* blob size*/

class Rstokeslet2D {
public:
	Rstokeslet2D(int numberOfParticles);
	virtual ~Rstokeslet2D();

	virtual void mkMatrix();
	virtual bool solve();
	virtual void getVelocities(int, double *iLoc, double *iVel);

	virtual double *getLocationArray();
	virtual double *getVelocityArray();
	virtual double *getForceArray();
	virtual double *getMatrix();

protected:
	double term1(double, double);
	double term2(double, double);

	int numberOfParticles; // number of particles

	double *loc;/* particle location */
	double *vel;/* particle velocity */
	double *foc;/* particle force */
	double *mat;/* 2Nx2N dense matrix */

private:	
	static const double tol = 1.0e-8;// tolarence for GS
	
	void matZero(); // zeros matrix
	void setElmMat(int i,int j,double elm); // insert a value at (i,j)
	void addElmMat(int i,int j,double elm); // add a value at (i,j)
};
