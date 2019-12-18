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


// main
/*int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << argv[0] << " [number of points]\n";
        return 0;
*/

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

    }
    
//    int num = std::atoi(argv[1]);
//    std::cout << "Number of Points=" << num << std::endl;
    
    double *phi = new double[num];
 //double *phi0 = new double[num];
 //   double *rhs = new double[num];

for (int n = 0; n < num; n++) {
    phi[n] = 0.0;
}
double dx = 1.0 / (num - 1);




// declare Sparse Matrix
//int size = num - 2; // size of inner points
//Eigen::SparseMatrix<double> matA(size * size);
//matA.reserve(Eigen::VectorXi::Constant(size * size, 5));
// vector for RHS
Eigen::VectorXd rhs(num);
// vector for solution
Eigen::VectorXd phi0(num);

// set up Matrix
double start = clock();

for (int i = 0; i < num; i++) {
    rhs[i] = -dx * dx / d * (dx * i);
    phi[i] = phi0[i] = 0.0;
}
// Boundary Conditions
phi[0] = phi0[0] = a;
phi[num - 1] = phi0[num - 1] = b;






double tcost = (clock() - start) / CLOCKS_PER_SEC;
std::cout << "Time cost setup = " << tcost << "(sec)\n";

// solve linear system with CG (Conjugate Gradient Method)
const double tol = 1.0e-8;
Eigen::ConjugateGradient < Eigen::SparseMatrix<double> > solver;
solver.setTolerance(tol);
solver.compute(phi);
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
