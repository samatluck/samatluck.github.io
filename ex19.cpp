//
// Large Scale Computing
// 2D Heat/Mass Transfer
// ex19.cpp :
// Solve for
// d2c/dx2 + d2c/dy2 + d2c/dz2 + 1 = 0
// With the boundary conditions of c = 0
// along lines of x=1,-1,  y=1,-1, and z=1,-1
//
// Parallel Jacobi

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <omp.h>

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << argv[0] << " [number of points]\n";
		return 0;
	}

	int num = std::atoi(argv[1]);
	std::cout << "Number of Points=" << num << std::endl;

	double *phi = new double[num * num * num];
	double *phi0 = new double[num * num * num];

	double *rhs = new double[num * num * num];
	double *ap = new double[num * num * num];
	double *ae = new double[num * num * num];
	double *aw = new double[num * num * num];
	double *an = new double[num * num * num];
	double *as = new double[num * num * num];
	double *at = new double[num * num * num];
	double *ab = new double[num * num * num];

	int *ie = new int[num * num * num];
	int *iw = new int[num * num * num];
	int *in = new int[num * num * num];
	int *is = new int[num * num * num];
	int *it = new int[num * num * num];
	int *ib = new int[num * num * num];

	/*assuming dx = dy : domain is 2x2 */
	double dx = 2.0 / (num - 1);

	/* Initialize */
#pragma omp parallel for schedule(runtime)
	for (int n = 0; n < num * num * num; n++) {
		phi[n] = 0.0;
		phi0[n] = 0.0;
	}

#pragma omp parallel for schedule(runtime)
	for (int n = 0; n < num * num * num; n++) {
		rhs[n] = -dx * dx;
		ap[n] = -6;
		ae[n] = 1;
		aw[n] = 1;
		an[n] = 1;
		as[n] = 1;
		at[n] = 1;
		ab[n] = 1;
		ie[n] = n + 1;
		iw[n] = n - 1;
		in[n] = n + num;
		is[n] = n - num;
		it[n] = n + num * num;
		ib[n] = n - num * num;
		if (((n % num) == 0) || ((n % num) == (num - 1)) || ((n / num) % num == 0)
				|| ((n / num) % num == (num - 1)) || ((n / (num * num)) == 0)
				|| ((n / (num * num)) == (num - 1))) {
			rhs[n] = phi[n] * ap[n];
			ae[n] = aw[n] = an[n] = as[n] = at[n] = ab[n] = 0;
			ie[n] = iw[n] = in[n] = is[n] = it[n] = ib[n] = n;
		}
	}

	/* computing for phi with Jacobi Method */
	int itc = 0;
	const double tol = 1.0e-8;
	double start = clock();
	double omp_start = omp_get_wtime();
	while (1) {
#pragma omp parallel for schedule(runtime)
		for (int n = 0; n < num * num  * num ; n++) {
			phi[n] = (rhs[n] - ae[n] * phi0[ie[n]] - aw[n] * phi0[iw[n]]
					- an[n] * phi0[in[n]] - as[n] * phi0[is[n]] - at[n] * phi0[it[n]] - ab[n] * phi0[ib[n]]) / ap[n];
		}
		double err = 0.0;
#pragma omp parallel for schedule(runtime) reduction(+:err)
		for (int n = 0; n < num * num  * num; n++) {
			double r = rhs[n] - ae[n] * phi[ie[n]] - aw[n] * phi[iw[n]]
					- an[n] * phi[in[n]] - as[n] * phi[is[n]] - at[n] * phi[it[n]] - ab[n] * phi[ib[n]] - ap[n] * phi[n];
			err += r * r;
			phi0[n] = phi[n];
		}
		itc++;
		if (err < tol * tol)
			break;
		if (itc % 1000 == 0) {
			std::cout << "# of Iteration=" << itc << " err=" << std::sqrt(err)
					<< std::endl;
		}
	}
	double tcost = (clock() - start) / CLOCKS_PER_SEC;
	double omp_tcost = omp_get_wtime() - omp_start;
	std::cout << "Number of Iteration=" << itc << std::endl;
	std::cout << "Time cost (CPU) = " << tcost << "(sec)\n";
	std::cout << "Time cost (WALL CLOCK)= " << omp_tcost << "(sec)\n";

	// Output Result
	std::ofstream ofile;
	ofile.open("ex19.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			for (int k = 0; k < num; k++) {
				double x = -1.0 + 2.0 * i / (num - 1);
				double y = -1.0 + 2.0 * j / (num - 1);
				double z = -1.0 + 2.0 * k / (num - 1);
				ofile << x << " " << y << " " << z << " "<< phi[i + num * j + num * num * k] << std::endl;
			}
		}
	}
	ofile.close();

	std::cout << "Done\n";
	return 0;
	printf("Done\n");
}

