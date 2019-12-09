//
//  Large Scale Computing
//  Stokes Flow in a Cavity
//  ex74.cpp
//
//  Use Lapack and Blas Libraries
//
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <new>
#include <sys/time.h>
#include "stokeslet2dLapack.h"

#define INTGRID 51    /* # of x-grid lines for internal velocity */

// timing method
double tsecond() {
	struct timeval tm;
	double t;
	static int base_sec = 0, base_usec = 0;

	gettimeofday(&tm, NULL);
	if (base_sec == 0 && base_usec == 0) {
		base_sec = tm.tv_sec;
		base_usec = tm.tv_usec;
		t = 0.0;
	} else {
		t = (double) (tm.tv_sec - base_sec) + ((double) (tm.tv_usec - base_usec)) / 1.0e6;
	}
	return t;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [Depth of Cavity]\n";
		return -1;
	}

	/* get inputed depth */
	double dp = std::atof(argv[1]);

	/* EPSILON = 0.005 defined by 'stokeslet2d.h' */	
	int numpdepth = (dp / EPSILON + 0.5);  // # of particles in depth 
	int numpwidth = (1.0 / EPSILON + 0.5); // # of particles in width 
	int numparticles = numpdepth * 2 + numpwidth * 2; // total # od particles 
	std::cout << "Total # of Particles= " << numparticles << std::endl;

	// Allocate Space 
	// Create an insetance of object Rstokeslet2D
	Rstokeslet2D *slet = new Rstokeslet2DLapack(numparticles);

	// get each arrays
	double *loc = slet->getLocationArray();
	double *vel = slet->getVelocityArray();

	/* set location & velocity of particles (blob)*/
	for (int i = 0; i < numparticles; i++) {
		if ((i >= 0) && (i <= numpwidth)) { /* top */
			loc[i * DIM] = -0.5 + EPSILON * i; /* x */
			loc[i * DIM + 1] = 0.0; /* y */
			vel[i * DIM] = 1.0;
			vel[i * DIM + 1] = 0.0;
		} else {
			if (i <= (numpwidth + numpdepth)) { /* right wall */
				loc[i * DIM] = 0.5; /* x */
				loc[i * DIM + 1] = -EPSILON * (i - numpwidth); /* y */
				vel[i * DIM] = 0.0;
				vel[i * DIM + 1] = 0.0;
			} else {
				if (i <= (2 * numpwidth + numpdepth)) { /* bottom */
					loc[i * DIM] = 0.5
							- EPSILON * (i - (numpwidth + numpdepth)); /* x */
					loc[i * DIM + 1] = -EPSILON * numpdepth; /* y */
					vel[i * DIM] = 0.0;
					vel[i * DIM + 1] = 0.0;
				} else { /* left wall */
					loc[i * DIM] = -0.5; /* x */
					loc[i * DIM + 1] = -EPSILON
							* ((2 * numpwidth + 2 * numpdepth) - i); /* y */
					vel[i * DIM] = 0.0;
					vel[i * DIM + 1] = 0.0;
				}
			}
		}
	}

	/* make 2Nx2N Matrix */
	double start = tsecond();
	slet->mkMatrix();
	double tcost = (tsecond() - start);
	std::cout << "Time cost for setup Matrix (wall-time) = " << tcost << "(sec)\n";

	/* Sovle linear ststem */
	start = tsecond();
	slet->solve();
	tcost = (tsecond() - start);
	std::cout << "Time cost to solve linear system (wall-time) = " << tcost << "(sec)\n";

	/* 
	 compute internal velocity 
	 */
	int nyg = (EPSILON * (numpdepth * (INTGRID - 1)));
	double *cvel = new double[INTGRID * nyg * DIM];
	double *cloc = new double[INTGRID * nyg * DIM];

	/* setting grid */
	for (int j = 0; j < nyg; j++) {
		for (int i = 0; i < INTGRID; i++) {
			cloc[DIM * (i + INTGRID * j)] = -0.5
					+ static_cast<double>(i) / (INTGRID - 1);
			cloc[DIM * (i + INTGRID * j) + 1] = -static_cast<double>(j)
					/ (INTGRID - 1);
		}
	}

	/* compute velocities */
	start = tsecond();
	slet->getVelocities(INTGRID * nyg, cloc, cvel);
	tcost = (tsecond() - start);
	std::cout << "Time cost to compute interior velocity (wall-time) = " << tcost
			<< "(sec)\n";

	/* out velocities */
	std::ofstream ofile;
	ofile.open("resEx74.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
	for (int j = 0; j < nyg; j++) {
		for (int i = 0; i < INTGRID; i++) {
			ofile << cloc[DIM * (i + INTGRID * j)] << " "
					<< cloc[DIM * (i + INTGRID * j) + 1] << " "
					<< cvel[DIM * (i + INTGRID * j)] << " "
					<< cvel[DIM * (i + INTGRID * j) + 1] << std::endl;
		}
	}
	ofile.close();

	/* free */
	delete[] cloc;
	delete[] cvel;
	delete slet;
}

