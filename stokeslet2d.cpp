/*
 Large Scale Computing
 Stokeslet for 2D
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "stokeslet2d.h"

//Constructor
Rstokeslet2D::Rstokeslet2D(int numberOfParticles) :
		numberOfParticles(numberOfParticles) {
	try {
		const int np = numberOfParticles;
		loc = new double[np * DIM];/* particle location */
		vel = new double[np * DIM];/* particle velocity */
		foc = new double[np * DIM];/* particle force */
		mat = new double[np * DIM * np * DIM];/* 2Nx2N dense matrix */
	} catch (std::bad_alloc& ba) {
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
		std::abort();
	}
}

//Destructor
Rstokeslet2D::~Rstokeslet2D() {
	delete[] mat;
	delete[] loc;
	delete[] vel;
	delete[] foc;
}

/* make Matrix */
void Rstokeslet2D::mkMatrix() {
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

/* Sovle Liear ststem with GS*/
bool Rstokeslet2D::solve() {
	const int np = numberOfParticles;
	int info;

	/*solve Ax=b */
	/* A (mat) is n x n matrix */
	int n = DIM * np;
	int itr = 0;
	while(1){
		// GS iteration step
		for (int i = 0 ; i < n ; i++){
			double diag = mat[i * n + i];
			double elm = vel[i];
			for (int j = 0 ; j < n ; j++){
				if (i != j){
					elm -=  mat[i * n + j] * foc[j];
				}
			}
			foc[i] = elm / diag;
		}
		
		// Check Convergence
		double err = 0.0;
		for (int i = 0 ; i < n ; i++){
			double elm = vel[i];
			for (int j = 0 ; j < n ; j++){
				elm -=  mat[i * n + j] * foc[j];
			}
			err += elm * elm;
		}		
		if (err < tol * tol) break;
		
		itr++;
		if (itr % 1000 == 0){
			std::cout << "R=" << std::sqrt(err) << std::endl;
		}
	}
	
	return true;
}

/* compute velocity */
void Rstokeslet2D::getVelocities(int nump, double *cloc, double *cvel) {
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

double *Rstokeslet2D::getLocationArray() {
	return loc;
}

double *Rstokeslet2D::getVelocityArray() {
	return vel;
}

double *Rstokeslet2D::getForceArray() {
	return foc;
}

double *Rstokeslet2D::getMatrix() {
	return mat;
}

double Rstokeslet2D::term1(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

double Rstokeslet2D::term2(double r, double ep) {
	double sq;

	sq = sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}

// zeros matrix
void Rstokeslet2D::matZero(){
	const int si = DIM * DIM * numberOfParticles * numberOfParticles;
	for (int i = 0; i < si; i++){
		mat[i] = 0.0;
	}
}

// insert a value at (i,j)
void Rstokeslet2D::setElmMat(int i,int j,double elm){
	const int offset = DIM * numberOfParticles;
	mat[i + j * offset] = elm;
}

// add a value at (i,j)
void Rstokeslet2D::addElmMat(int i,int j,double elm){
	const int offset = DIM * numberOfParticles;
	mat[i + j * offset] += elm;
}