//
//  Large Scale Computing
//  Heat/Mass Transfer
//  ex01.cpp
// Jacobi Method
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
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

    
    
        }
    }
double tcost = (clock() - start) / CLOCKS_PER_SEC;
std::cout << "Time cost setup = " << tcost << "(sec)\n";

const double tol = 1.0e-8;
Eigen::ConjugateGradient < Eigen::SparseMatrix<double> > solver;
solver.setTolerance(tol);
solver.compute(phi0);
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
