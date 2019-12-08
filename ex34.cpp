//
// Large Scale Computing
// 2D Heat/Mass Transfer
// ex34.cpp 
// Solve for
// d2c/dx2 + d2c/dy2 + 1 = 0
// With the boundary conditions of c = 0 
// along lines of x=1,-1, and y=1, and -1
//
//  Conjugate Gradient Method with Eigen3
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// index for all points
int gIdx(int i, int j, int num) {
	return i + num * j;
}

// index for interior points
int iIdx(int i, int j, int num) {
	return i - 1 + (num - 2) * (j - 1);
}

// main
int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << argv[0] << " [number of points]\n";
		return 0;
	}

	int num = std::atoi(argv[1]);
	std::cout << "Number of Points=" << num << std::endl;

	double *phi = new double[num * num];
	// Initialize
	for (int n = 0; n < num * num; n++) {
		phi[n] = 0.0;
	}

	// assuming dx = dy : domain is 2x2
	double dx = 2.0 / (num - 1);

	// declare Sparse Matrix
	int size = num - 2; // size of inner points
	Eigen::SparseMatrix<double> matA(size * size, size * size);
	matA.reserve(Eigen::VectorXi::Constant(size * size, 5));
	// vector for RHS
	Eigen::VectorXd rhs(size * size);
	// vector for solution
	Eigen::VectorXd sol(size * size);

	// set up Matrix
	double start = clock();
	for (int i = 1; i < num - 1; i++) {
		for (int j = 1; j < num - 1; j++) {

			int n = iIdx(i, j, num);

			// initial guess
			sol(n) = phi[gIdx(i, j, num)];

			// RHS
			rhs(n) = -dx * dx;

			// diagonal element
			matA.coeffRef(n, n) = -4.0;

			// east
			if (i != num - 2) {
				matA.coeffRef(n, iIdx(i + 1, j, num)) = 1.0;
			} else {
				rhs(n) -= phi[gIdx(i + 1, j, num)];
			}

			// west
			if (i != 1) {
				matA.coeffRef(n, iIdx(i - 1, j, num)) = 1.0;
			} else {
				rhs(n) -= phi[gIdx(i - 1, j, num)];
			}

			// north
			if (j != num - 2) {
				matA.coeffRef(n, iIdx(i, j + 1, num)) = 1.0;
			} else {
				rhs(n) -= phi[gIdx(i, j + 1, num)];
			}

			// south
			if (j != 1) {
				matA.coeffRef(n, iIdx(i, j - 1, num)) = 1.0;
			} else {
				rhs(n) -= phi[gIdx(i, j - 1, num)];
			}
		}
	}
	double tcost = (clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Time cost setup = " << tcost << "(sec)\n";

	// solve linear system with CG (Conjugate Gradient Method)
	const double tol = 1.0e-8;
	Eigen::ConjugateGradient < Eigen::SparseMatrix<double> > solver;
	solver.setTolerance(tol);
	solver.compute(matA);
	sol = solver.solve(rhs);

	double tcost2 = (clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Time cost solver = " << tcost2 - tcost << "(sec)\n";
	std::cout << "#iterations:     " << solver.iterations() << std::endl;
	std::cout << "estimated error: " << solver.error() << std::endl;

	// restore solution
	for (int i = 1; i < num - 1; i++) {
		for (int j = 1; j < num - 1; j++) {
			phi[gIdx(i, j, num)] = sol(iIdx(i, j, num));
		}
	}

	// Output Result
	std::ofstream ofile;
	ofile.open("res2d.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
	for (int i = 0; i < num; i++) {
		for (int j = 0; j < num; j++) {
			double x = -1.0 + 2.0 * i / (num - 1);
			double y = -1.0 + 2.0 * j / (num - 1);
			ofile << x << " " << y << " " << phi[gIdx(i, j, num)] << std::endl;
		}
	}
	ofile.close();

	std::cout << "Done\n";
	return 0;
}
