/*
 Large Scale Computing
 Stokeslet for 2D
 MPI Version
 
 Ricardo Cortez,"The Method of Regularized Stokeslets", 
 (2001), SIAM J. Sci. Comput., Vol.23, No.4, pp.1204
 
 assuming mu = 1
 
 Use MPI
*/
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include "stokeslet2d_dist.h"
//Constructor
Rstokeslet2D_dist::Rstokeslet2D_dist(int numberOfParticles) :
		numberOfParticles(numberOfParticles) {

	numproc = MPI::COMM_WORLD.Get_size();
	myid = MPI::COMM_WORLD.Get_rank();

	/* divide blobs into # of cpus */
	myStart = (numberOfParticles / numproc) * myid;
	if (numberOfParticles % numproc > myid) {
		myStart += myid;
		myEnd = myStart + (numberOfParticles / numproc) + 1;
	} else {
		myStart += numberOfParticles % numproc;
		myEnd = myStart + (numberOfParticles / numproc);
	}
	mySize = myEnd - myStart;

	/* collect start end info */
	oStart = new int[numproc];
	oEnd = new int[numproc];
	MPI_Allgather(&mySize, 1, MPI_INT, oEnd, 1, MPI_INT, MPI_COMM_WORLD);
	oStart[0] = 0;
	for (int cpu = 1; cpu < numproc; cpu++) {
		oStart[cpu] = oStart[cpu - 1] + oEnd[cpu - 1];
		oEnd[cpu - 1] += oStart[cpu - 1];
	}
	oEnd[numproc - 1] += oStart[numproc - 1];

	try {
		const int np = numberOfParticles;
		loc = new double[mySize * DIM];/* particle location */
		vel = new double[mySize * DIM];/* particle velocity */
		foc = new double[mySize * DIM];/* particle force */
		mat = new double[mySize * DIM * numberOfParticles * DIM];/* 2Nx2N dense matrix */
	} catch (std::bad_alloc& ba) {
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
		std::abort();
	}
}

//Destructor
Rstokeslet2D_dist::~Rstokeslet2D_dist() {
	delete[] mat;
	delete[] loc;
	delete[] vel;
	delete[] foc;
}

/* make Matrix */
void Rstokeslet2D_dist::mkMatrix() {
	const double ep = EPSILON;
	const int np = numberOfParticles;
	
	/* get max size */
	int maxsize = 0;
	for (int cpu = 0; cpu < numproc; cpu++) {
		maxsize = std::max(oEnd[cpu] - oStart[cpu], maxsize);
	}

	/* alloc */
	double *oloc = new double[maxsize * DIM];

	/*zeros mat*/
	for (int i = 0; i < DIM * DIM * np * mySize; i++) {
		mat[i] = 0.0;
	}

	/* make mat */
	for (int cpu = 0; cpu < numproc; cpu++) {
		/* get location vector from others */
		double *ploc = oloc;
		if (myid == cpu) {
			ploc = loc;
		}
		MPI_Bcast(ploc, DIM * (oEnd[cpu] - oStart[cpu]), MPI_DOUBLE, cpu,
				MPI_COMM_WORLD);

		/* compute local matrix*/
		for (int i = oStart[cpu]; i < oEnd[cpu]; i++) {
			for (int j = oStart[myid]; j < oEnd[myid]; j++) {
				int li = i - oStart[cpu];
				int lj = j - oStart[myid];
				double dx = ploc[li * DIM] - loc[lj * DIM];
				double dy = ploc[li * DIM + 1] - loc[lj * DIM + 1];
				double r = std::sqrt(dx * dx + dy * dy);

				double tr1 = term1(r, ep) / (4.0 * M_PI);
				double tr2 = term2(r, ep) / (4.0 * M_PI);

				mat[i * DIM + lj * DIM * np * DIM] += -tr1 + tr2 * dx * dx;
				mat[i * DIM + (lj * DIM + 1) * np * DIM] += tr2 * dx * dy;
				mat[i * DIM + 1 + lj * DIM * np * DIM] += tr2 * dx * dy;
				mat[i * DIM + 1 + (lj * DIM + 1) * np * DIM] += -tr1
						+ tr2 * dy * dy;
			}
		}
	}

	delete[] oloc;
}

/* Sovle Liear ststem */
/* CG */
bool Rstokeslet2D_dist::solve() {

	int n = DIM * numberOfParticles;
	double *b = vel;
	double *x = foc;
	int ms = myStart * DIM;
	int me = myEnd * DIM;

	/* Allocate Working Spaces */
	double *r = new double[mySize * DIM];
	double *p = new double[mySize * DIM];
	double *Ap = new double[mySize * DIM];
	double *work = new double[n];

	/* initialize */
	for (int i = 0; i < mySize * DIM; i++) {
		r[i] = p[i] = Ap[i] = 0.0;
	}

	/* |b| */
	double w = 0.0;
	for (int i = 0; i < mySize * DIM; i++) {
		w += b[i] * b[i]; /* w = b.b */
	}
	double normb;
	MPI_Allreduce(&w, &normb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	normb = std::sqrt(normb);

	/* compute r0 = b - Ax */
	for (int i = 0; i < n; i++)
		work[i] = 0.0;
	for (int i = ms; i < me; i++) {
		int li = i - ms;
		work[i] = b[li];
	}
	for (int j = ms; j < me; j++) {
		int lj = j - ms;
		for (int i = 0; i < n; i++) {
			work[i] -= mat[i + lj * n] * x[lj];
		}
	}
	for (int cpu = 0; cpu < numproc; cpu++) {
		MPI_Reduce(&work[DIM * oStart[cpu]], r, DIM * (oEnd[cpu] - oStart[cpu]),
				MPI_DOUBLE, MPI_SUM, cpu, MPI_COMM_WORLD);
	}

	w = 0.0;
	for (int i = 0; i < me - ms; i++) {
		p[i] = r[i]; /* p = r */
		w += r[i] * r[i]; /* rr = r.r */
	}
	double rr;
	MPI_Allreduce(&w, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	/* cg iteration */
	int count = 0;
	while (rr > tol * tol * normb * normb) {
		/* Ap = A*p */
		for (int i = 0; i < n; i++)
			work[i] = 0.0;
		for (int j = ms; j < me; j++) {
			int lj = j - ms;
			for (int i = 0; i < n; i++) {
				work[i] += mat[i + lj * n] * p[lj];
			}
		}
		for (int cpu = 0; cpu < numproc; cpu++) {
			MPI_Reduce(&work[DIM * oStart[cpu]], Ap,
					DIM * (oEnd[cpu] - oStart[cpu]), MPI_DOUBLE, MPI_SUM, cpu,
					MPI_COMM_WORLD);
		}

		/* pAp = p.Ap */
		w = 0.0;
		for (int i = 0; i < me - ms; i++) {
			w += p[i] * Ap[i];
		}
		double pAp;
		MPI_Allreduce(&w, &pAp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		/* alpha = r.r / p.Ap */
		double alpha = rr / pAp;

		/* Beta */
		for (int i = 0; i < me - ms; i++) {
			x[i] += alpha * p[i]; /* x += alpha * p */
			r[i] -= alpha * Ap[i]; /* r -= alpha * Ap */
		}

		/* rr1 = r.r */
		w = 0.0;
		for (int i = 0; i < me - ms; i++) {
			w += r[i] * r[i];
		}
		double rr1;
		MPI_Allreduce(&w, &rr1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		double beta = rr1 / rr;

		/* p = r + beta * p */
		for (int i = 0; i < me - ms; i++) {
			p[i] = r[i] + beta * p[i];
		}

		rr = rr1;
		count++;
	}

	/* Deallocate Working Spaces */
	delete[] r;
	delete[] p;
	delete[] Ap;
	delete[] work;
	
	return true;
}

/* compute velocity */
void Rstokeslet2D_dist::getVelocities(int nump, double *cloc, double *cvel) {
	const double ep = EPSILON;
	const int np = numberOfParticles;
	double ov[2];

	for (int j = 0; j < nump; j++) {

		/* zeros */
		ov[0] = 0.0;
		ov[1] = 0.0;

		/* loop for partilces (blob) */
		for (int i = myStart; i < myEnd; i++) {
			int li = i - myStart;
			double dx = cloc[j * DIM] - loc[li * DIM];
			double dy = cloc[j * DIM + 1] - loc[li * DIM + 1];
			double r = sqrt(dx * dx + dy * dy);

			double tr1 = term1(r, ep) / (4.0 * M_PI);
			double tr2 = term2(r, ep) / (4.0 * M_PI);

			tr2 *= foc[li * DIM] * dx + foc[li * DIM + 1] * dy;

			ov[0] += -foc[li * DIM] * tr1 + tr2 * dx;
			ov[1] += -foc[li * DIM + 1] * tr1 + tr2 * dy;
		}

		/* sum up to CPU0 */
		MPI_Reduce(ov, &cvel[j * DIM], 2, MPI_DOUBLE, MPI_SUM, 0,
				MPI_COMM_WORLD);
	}
}

double *Rstokeslet2D_dist::getLocationArray() {
	return loc;
}

double *Rstokeslet2D_dist::getVelocityArray() {
	return vel;
}

double *Rstokeslet2D_dist::getForceArray() {
	return foc;
}

double *Rstokeslet2D_dist::getMatrix() {
	return mat;
}

double Rstokeslet2D_dist::term1(double r, double ep) {
	double sq;

	sq = std::sqrt(r * r + ep * ep);

	return std::log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

double Rstokeslet2D_dist::term2(double r, double ep) {
	double sq;

	sq = std::sqrt(r * r + ep * ep);

	return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}
