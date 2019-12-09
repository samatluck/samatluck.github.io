/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use Eigen3
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "stokeslet2dEigen3.h"

//Constructor
Rstokeslet2DEigen3::Rstokeslet2DEigen3(int np) :
		numberOfParticles(np), loc(np * DIM), vel(np * DIM), foc(np * DIM), mat(
				np * DIM, np * DIM) {
}

//Destructor
Rstokeslet2DEigen3::~Rstokeslet2DEigen3() {

}
/* make Matrix */
void Rstokeslet2DEigen3::mkMatrix() {
	const int np = numberOfParticles;
	const double ep = EPSILON;

	/*zeros mat*/
	matZero();

	/* make mat */
	/* Because of symmetric, we can loop over half */
	for (int i = 0; i < np; i++) {
		for (int j = i; j < np; j++) {
			double dx = loc(i * DIM) - loc(j * DIM);
			double dy = loc(i * DIM + 1) - loc(j * DIM + 1);
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
bool Rstokeslet2DEigen3::solve() {
	
	Eigen::PartialPivLU < Eigen::MatrixXd > solver(mat);
	foc = solver.solve(vel);

	return true;
}

/* compute velocity */
void Rstokeslet2DEigen3::getVelocities(int nump, double *cloc, double *cvel) {
	const int np = numberOfParticles;
	const double ep = EPSILON;

	/* loop for grid */
	for (int p = 0; p < nump; p++) {
		/* zeros */
		cvel[p * DIM] = 0.0;
		cvel[p * DIM + 1] = 0.0;

		/* loop for partilces (blob) */
		for (int i = 0; i < np; i++) {
			double dx = cloc[p * DIM] - loc(i * DIM);
			double dy = cloc[p * DIM + 1] - loc(i * DIM + 1);
			double r = sqrt(dx * dx + dy * dy);

			double tr1 = term1(r, ep) / (4.0 * M_PI);
			double tr2 = term2(r, ep) / (4.0 * M_PI);

			tr2 *= foc(i * DIM) * dx + foc(i * DIM + 1) * dy;

			cvel[p * DIM] += -foc(i * DIM) * tr1 + tr2 * dx;
			cvel[p * DIM + 1] += -foc(i * DIM + 1) * tr1 + tr2 * dy;
		}
	}
}

double *Rstokeslet2DEigen3::getLocationArray() {
	return loc.data();
}

double *Rstokeslet2DEigen3::getVelocityArray() {
	return vel.data();
}

double *Rstokeslet2DEigen3::getForceArray() {
	return foc.data();
}

double *Rstokeslet2DEigen3::getMatrix() {
	return mat.data();
}

double Rstokeslet2DEigen3::term1(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

double Rstokeslet2DEigen3::term2(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}

// zeros matrix
void Rstokeslet2DEigen3::matZero() {
	mat.setZero();
}

// insert a value at (i,j)
void Rstokeslet2DEigen3::setElmMat(int i, int j, double elm) {
	mat(i, j) = elm;
}

// add a value at (i,j)
void Rstokeslet2DEigen3::addElmMat(int i, int j, double elm) {
	mat(i, j) += elm;
}