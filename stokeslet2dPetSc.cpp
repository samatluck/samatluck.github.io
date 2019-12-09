/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use PetSc
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "stokeslet2dPetSc.h"

//Constructor
Rstokeslet2DPetSc::Rstokeslet2DPetSc(int np) :
		numberOfParticles(np) {
	numproc = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();	
	
	int myStart = (np / numproc) * myid;
	int myend;
	if (np % numproc > myid) {
		myStart += myid;
		myEnd = myStart + (np / numproc) + 1;
	} else {
		myStart += np % numproc;
		myEnd = myStart + (np / numproc);
	}
	int mySize = myEnd - myStart;

	
	VecCreateMPI(PETSC_COMM_WORLD, mySize * DIM, np * DIM, &loc);
	VecCreateMPI(PETSC_COMM_WORLD, mySize * DIM, np * DIM, &vel);
	VecCreateMPI(PETSC_COMM_WORLD, mySize * DIM, np * DIM, &foc);

#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3)
	MatCreateAIJ(PETSC_COMM_WORLD, mySize * DIM, mySize * DIM, np * DIM, np * DIM, PETSC_DECIDE, NULL, PETSC_DECIDE, NULL, &mat);
#else
	MatCreateMPIAIJ(PETSC_COMM_WORLD, mySize * DIM, mySize * DIM, np * DIM, np * DIM, PETSC_DECIDE, NULL, PETSC_DECIDE, NULL, &mat);
#endif

}

//Destructor
Rstokeslet2DPetSc::~Rstokeslet2DPetSc() {

}
/* make Matrix */
void Rstokeslet2DPetSc::mkMatrix() {
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
bool Rstokeslet2DPetSc::solve() {
	
	Eigen::PartialPivLU < Eigen::MatrixXd > solver(mat);
	foc = solver.solve(vel);

	return true;
}

/* compute velocity */
void Rstokeslet2DPetSc::getVelocities(int nump, double *cloc, double *cvel) {
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

double *Rstokeslet2DPetSc::getLocationArray() {
	return loc.data();
}

double *Rstokeslet2DPetSc::getVelocityArray() {
	return vel.data();
}

double *Rstokeslet2DPetSc::getForceArray() {
	return foc.data();
}

double *Rstokeslet2DPetSc::getMatrix() {
	return mat.data();
}

double Rstokeslet2DPetSc::term1(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

double Rstokeslet2DPetSc::term2(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}

// zeros matrix
void Rstokeslet2DPetSc::matZero() {
	mat.setZero();
}

// insert a value at (i,j)
void Rstokeslet2DPetSc::setElmMat(int i, int j, double elm) {
	mat(i, j) = elm;
}

// add a value at (i,j)
void Rstokeslet2DPetSc::addElmMat(int i, int j, double elm) {
	mat(i, j) += elm;
}