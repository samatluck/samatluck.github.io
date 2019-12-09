//
//  Large Scale Computing
//  Heat/Mass Transfer
//  ex01.cpp
// Jacobi Method
//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main(int argc, char **argv) {
	int num;
	double a, b, d;
	//  set parameters
	std::cout << "Number of Grid Points [N]\n";
	std::cout << "Boundary Value at Left [A]\n";
	std::cout << "Boundary Value at Right [B]\n";
	std::cout << "Diffusion Coefficient [D]\n";
	std::cout << "Input [N] [A] [B] [D] : ";
	std::cin >> num >> a >> b >> d;

	if (std::cin.fail()) {
		std::cout << "BYE\n";
		return 0;
	}

	// Allocation Array
	double *phi = new double[num];
	double *phi0 = new double[num];
	double *rhs = new double[num];

	// Setup
	double dx = 1.0 / (num - 1);
	for (int i = 0; i < num; i++) {
		rhs[i] = -dx * dx / d * (dx * i);
		phi[i] = phi0[i] = 0.0;
	}

	// Boundary Conditions
	phi[0] = phi0[0] = a;
	phi[num - 1] = phi0[num - 1] = b;

	//
	// Solve with Jacobi Method
	int itc = 0;
	const double tol = 1.0e-8;
	while (1) {
		// update phi
		for (int i = 1; i < num - 1; i++) {
			phi[i] = -0.5 * (rhs[i] - phi0[i + 1] - phi0[i - 1]);
		}

		// Check Convergence
		double merr = 0.0;
		for (int i = 1; i < num - 1; i++) {
			double r = rhs[i] - phi[i + 1] - phi[i - 1] + 2.0 * phi[i];
			merr += r * r;
			phi0[i] = phi[i];
		}
		merr = std::sqrt(merr);

		itc++;
		if (merr < tol)
			break;

		if ((itc % 10000) == 0) {
			std::cout << "itc=" << itc << " err=" << merr << std::endl;
		}
	}
	std::cout << "Number of Iteration=" << itc << std::endl;

	// Output Result
	std::ofstream ofile;
	ofile.open("ex01.dat");
	ofile << std::setprecision(16);
	ofile << std::scientific;
	for (int i = 0; i < num; i++) {
		ofile << dx * i << " " << phi[i] << std::endl;
	}
	ofile.close();

	std::cout << "Done\n";
	return 0;
}
