/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 */
#define DIM 2
#define EPSILON 0.005 /* blob size*/

class Rstokeslet2D_dist {
public:
	Rstokeslet2D_dist(int numberOfParticles);
	virtual ~Rstokeslet2D_dist();

	void mkMatrix();
	bool solve();
	void getVelocities(int nump, double *iLoc, double *iVel);

	double *getLocationArray();
	double *getVelocityArray();
	double *getForceArray();
	double *getMatrix();

	int myStart;
	int myEnd;
	int mySize;
	
private:
	double term1(double, double);
	double term2(double, double);

	int numberOfParticles; // number of particles

	double *loc;/* particle location */
	double *vel;/* particle velocity */
	double *foc;/* particle force */
	double *mat;/* 2Nx2N dense matrix */

	int *oStart;
	int *oEnd;

	int myid;
	int numproc;
		
	static const double tol = 1.0e-8; // tolarence for CG
};

