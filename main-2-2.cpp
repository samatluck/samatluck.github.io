#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

// index for all points
int (int i, int num) {
    return i + num;
}

// index for interior points
int (int i, int num) {
    return i - 1 + (num - 2);
}

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
    double *phi = new double[num];
    double *phi0 = new double[num];
    
    double dx = 1.0 / (num - 1);
    for (int i = 0; i < num; i++) {
        phi[i] = phi0[i] = ((-(dx*dx)/d)*exp(i));
    }
    phi[0] = phi0[0] = a;
    phi[num - 1] = phi0[num - 1] = b;
    
    // declare Sparse Matrix
    int size = num - 2; // size of inner points
    Eigen::SparseMatrix<double> matA(size, size);
    matA.reserve(Eigen::VectorXi::Constant(size, 3));
    // vector for RHS
    Eigen::VectorXd rhs(size);
    // vector for solution
    Eigen::VectorXd sol(size);
    
    // set up Matrix
    double start = clock();
    for (int i = 1; i < num - 1; i++) {
        int n= i-1;
        sol(n) = phi[(i, num)];
        
        // RHS
        rhs(n) = (-dx * dx) ;
        
        // diagonal element
        matA.coeffRef(n,n) = -2.0;
        
        // east
        if (i != num - 2) {
            matA.coeffRef(n, n+1) = 1.0;
        } else {
            rhs(n) -= ((-(dx*dx)/d)*exp(i)) - b;
        }
        
        // west
        if (i != 1) {
            matA.coeffRef(n, n-1) = 1.0;
        } else {
            rhs(n) -= ((-(dx*dx)/d)*exp(i)) - a;
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
    for (int i = 0; i < num - 1; i++) {
        phi[(i, num)] = sol((i, num));
    }
    
    
    // Output Result
    std::ofstream ofile;
    ofile.open("res2d.dat");
    ofile << std::setprecision(16);
    ofile << std::scientific;
    for (int i = 1; i < num; i++) {
        double x = -1.0 + 2.0 / (num - 1);
        double y = -1.0 + 2.0 / (num - 1);
        ofile << x << " " << y << " " << phi[(i, num)] << std::endl;
    }
    ofile.close();
    
    std::cout << "Done\n";
    return 0;
}
