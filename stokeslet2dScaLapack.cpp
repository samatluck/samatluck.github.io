/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use Scalapack
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include "stokeslet2dScaLapack.h"

//Constructor
Rstokeslet2DScalapack::Rstokeslet2DScalapack(int np, int nump, int numq) :
		numberOfParticles(np), sclwrap(np * DIM, BLOCK_SIZE, nump, numq) {
	try {
		const int np = numberOfParticles;
		loc = new double[np * DIM];/* particle location */
		vel = new double[np * DIM];/* particle velocity */
		foc = new double[np * DIM];/* particle force */
	} catch (std::bad_alloc& ba) {
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
		std::abort();
	}
	numproc = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();	
}

//Destructor
Rstokeslet2DScalapack::~Rstokeslet2DScalapack() {

}
/* make Matrix */
void Rstokeslet2DScalapack::mkMatrix() {
	const int np = numberOfParticles;
	const double ep = EPSILON;

	/*zeros mat*/
	matZero();

	/* make mat */
	/* Because of symmetric, we can loop over half */
	for (int i = 0; i < np; i++) {
		for (int j = i; j < np; j++) {
			double dx = loc[i * DIM] - loc[j * DIM];
			double dy = loc[i * DIM + 1] - loc[j * DIM + 1];
			double r = std::sqrt(dx * dx + dy * dy);

			double tr1 = term1(r, ep) / (4.0 * M_PI);
			double tr2 = term2(r, ep) / (4.0 * M_PI);

			/* diagoanl elements and upper half */
			double elm11 = -tr1 + tr2 * dx * dx;
			double elm12 = tr2 * dx * dy;
			double elm21 = tr2 * dx * dy;
			double elm22 = -tr1 + tr2 * dy * dy;
			addElmMat(i * DIM, j * DIM, elm11);
			addElmMat(i * DIM, j * DIM + 1, elm12);
			addElmMat(i * DIM + 1, j * DIM, elm21);
			addElmMat(i * DIM + 1, j * DIM + 1, elm22);
			if (i != j) {
				/* lower half */
				addElmMat(j * DIM, i * DIM, elm11);
				addElmMat(j * DIM, i * DIM + 1, elm12);
				addElmMat(j * DIM + 1, i * DIM, elm21);
				addElmMat(j * DIM + 1, i * DIM + 1, elm22);
			}
		}
	}
}

/* Sovle Liear ststem with LU*/
bool Rstokeslet2DScalapack::solve() {
	const int si = DIM * numberOfParticles;
	for (int i = 0; i < si; i++) {
		sclwrap.setValueRHS(i, vel[i]);
	}
	sclwrap.pdgesv();

	for (int i = 0; i < si; i++) {
		foc[i] = sclwrap.getValueRHS(i);
	}

	return true;
}

void Rstokeslet2DScalapack::collectResults(){
	const int si = DIM * numberOfParticles;
	double *ffoc = new double[si];
	MPI::COMM_WORLD.Allreduce(foc, ffoc, si, MPI::DOUBLE, MPI::SUM);
	for (int i = 0 ; i < si ; i++){
		foc[i] = ffoc[i];
	}
	delete [] ffoc;
}

/* compute velocity */
void Rstokeslet2DScalapack::getVelocities(int nump, double *cloc,
		double *cvel) {
	const int np = numberOfParticles;
	const double ep = EPSILON;

	/* loop for grid */
	for (int p = 0; p < nump; p++) {
		/* zeros */
		cvel[p * DIM] = 0.0;
		cvel[p * DIM + 1] = 0.0;

		/* loop for partilces (blob) */
		for (int i = 0; i < np; i++) {
			double dx = cloc[p * DIM] - loc[i * DIM];
			double dy = cloc[p * DIM + 1] - loc[i * DIM + 1];
			double r = sqrt(dx * dx + dy * dy);

			double tr1 = term1(r, ep) / (4.0 * M_PI);
			double tr2 = term2(r, ep) / (4.0 * M_PI);

			tr2 *= foc[i * DIM] * dx + foc[i * DIM + 1] * dy;

			cvel[p * DIM] += -foc[i * DIM] * tr1 + tr2 * dx;
			cvel[p * DIM + 1] += -foc[i * DIM + 1] * tr1 + tr2 * dy;
		}
	}
}

double *Rstokeslet2DScalapack::getLocationArray() {
	return loc;
}

double *Rstokeslet2DScalapack::getVelocityArray() {
	return vel;
}

double *Rstokeslet2DScalapack::getForceArray() {
	return foc;
}

double Rstokeslet2DScalapack::term1(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

double Rstokeslet2DScalapack::term2(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}

// zeros matrix
void Rstokeslet2DScalapack::matZero() {
	const int si = DIM * numberOfParticles;
	for (int i = 0; i < si; i++) {
		for (int j = 0; j < si; j++) {
			sclwrap.setValue(i, j, 0.0);
		}
	}
}

// insert a value at (i,j)
void Rstokeslet2DScalapack::setElmMat(int i, int j, double elm) {
	sclwrap.setValue(i, j, elm);
}

// add a value at (i,j)
void Rstokeslet2DScalapack::addElmMat(int i, int j, double elm) {
	sclwrap.setValue(i, j, sclwrap.getValue(i, j) + elm);
}
