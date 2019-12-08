//
//  Large Scale Computing
//  Convective Heat/Mass Transfer
//  ex42.cpp : use Lapack 
//
//  duphi/dx = d2phi/dx2
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#ifdef MKL_ILP64
#include <mkl.h>
#include <mkl_cblas.h>
#include <mkl_lapack.h>
#else
extern "C" {
#include <cblas.h>
	int dgtsv_(int *, int *, double *, double *, double *, double *, int *, int *);
}
#endif
#endif

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
	if (argc < 6) {
		printf("Usage:%s [NUM] [A] [B] [D] [U]\n", argv[0]);
		return 0;
	}

	//set parameters
#ifdef MKL_ILP64
	MKL_INT num = std::atoi(argv[1]);
#else
	int num = std::atoi(argv[1]);
#endif
	double a = std::atof(argv[2]);
	double b = std::atof(argv[3]);
	double d = std::atof(argv[4]);
	double u = std::atof(argv[5]);

	std::cout << "num=" << num << " A=" << a << " B=" << b << " D=" << d << " U=" << u << std::endl;

	// Memory Allocation
#ifdef MKL_ILP64
	double *phi = (double *)mkl_malloc(num * sizeof(double),64);
#else
	double *phi = new double[num];
#endif
	// Setup
	double dx = 1.0 / (num - 1);
	for (int i = 0; i < num; i++) {
		phi[i] = 0.0;
	}

	double ae = d / dx + std::max(-u, 0.0);
	double aw = d / dx + std::max(u, 0.0);
	double ap = ae + aw;

	/* making tri-digonal matrix */
#ifdef MKL_ILP64
	double *diag = (double *)mkl_malloc(num * sizeof(double),64);
	double *udiag = (double *)mkl_malloc((num - 1) * sizeof(double),64);
	double *ldiag = (double *)mkl_malloc((num - 1) * sizeof(double),64);
#else
	double *diag = new double[num];
	double *udiag = new double[num - 1];
	double *ldiag = new double[num - 1];
#endif
	for (int i = 0; i < num; i++) {
		diag[i] = -ap;
		if (i != num - 1) {
			udiag[i] = ae;
		}
		if (i != 0) {
			ldiag[i - 1] = aw;
		}
	}

	/* RHS & Boundary Condition*/
	for (int i = 0; i < num; i++) {
		phi[i] = 0.0; /* use phi for RHS array */
	}
	phi[0] -= a * aw;
	phi[num - 1] -= b * ae;

	/* call Lapack */
#ifdef MKL_ILP64
	MKL_INT nrhs = 1; /* # of RHS = 1*/
	MKL_INT info;
#else
	int nrhs = 1; /* # of RHS = 1*/
	int info;
#endif
	double start = tsecond();
	dgtsv_(&num, &nrhs, ldiag, diag, udiag, phi, &num, &info);
	double tcost = (tsecond() - start);
	std::cout << "Time cost (CPU) = " << tcost << "(sec)\n";

	if (info == 0) {
		std::cout << "successfully done\n";
	}
	if (info < 0) {
		std::cout << "the " << -info << "-th argument had an illegal value\n";
		return -1;
	}
	if (info > 0) {
		std::cout << "U(" << info << "," << info << ") is exactly zero.\n";
		return -1;
	}

	// Output Result
	std::ofstream ofile;
	ofile.open("resEx40.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
	for (int i = 0; i < num; i++) {
		ofile << dx * i << " " << phi[i] << std::endl;
	}
	ofile.close();

	delete[] phi;
	std::cout << "Done\n";
	return 0;
}
