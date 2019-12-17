#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <omp.h>


int main(int argc, const char * argv[]) {
    if (argc < 2) {
        std::cout << argv[0] << " [size]\n";
        return 0;
    }
#pragma omp declare target
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
        
    int size = std::atoi(argv[1]);
    std::cout << "Size=" << size << std::endl;
    
    // allocation
    double *xvec = new double[size];
    double *bvec = new double[size];
    // Make 2D array
    double *amatBlock = new double[size * size];
    double **amat = new double*[size];
    for (int i = 0 ; i < size ; i++){
        amat[i] = amatBlock + size * i;
    }
    
    // setup matrix A lower elements
    for (int i = 0; i < size; i++) {
        double diag = 0.0;
        for (int j = 0; j < size; j++) {
            double elm = 1 + (rand() % 100) * 0.01;
            amat[i][j] = elm;
            diag += elm;
        }
        // |A_ii| >= sum(A_ij) i != j
        amat[i][i] = -diag;
    }
    
    // setup vector x elements
    for (int i = 0; i < size; i++) {
        xvec[i] = i;
    }
    
    // compute rhs
    std::fill(bvec,bvec+size,0.0);
    for (int i = 0; i < size ; i++) {
        for (int j = 0; j < size ; j++) {
            bvec[i] += amat[i][j] * xvec[j];
        }
    }
    
    // compute sol = A^-1 b with JACOBI
    double *sol = new double[size]; // iteration step k+1
    double *sol0 = new double[size];// iteration step k
    std::fill(sol,sol+size,0.0);
    std::fill(sol0,sol0+size,0.0);
    const double tol = 1.0e-8;
    int count = 0;
    double start = clock();
    double omp_start = omp_get_wtime();
    while(1){
        // JACOBI step x(k+1) = A_ii^-1 (b - A_ij (i!=j) x(k))
        for (int i = 0 ; i < size ; i++){
            double sumAs = 0.0;
            for (int j = 0 ; j < size ; j++){
                if (j == i) continue;
                sumAs += amat[i][j] * sol0[j];
            }
            sol[i] = (bvec[i] - sumAs) / amat[i][i];
        }
        
        // compute residual |r| = |b-Ax|
        double r = 0.0;
#pragma omp end declare target
#pragma omp for schedule(dynamic,10) 
        for (int i = 0 ; i < size ; i++){
            double sumAs = bvec[i];
            for (int j = 0 ; j < size ; j++){
                sumAs -= amat[i][j] * sol[j];
            }
            r += sumAs * sumAs;
            sol0[i] = sol[i];
            
        }
        r = std::sqrt(r);
        count++;
        if (r < tol) break;
        
        if (count % 1000 == 0){
            std::cout << "count=" << count << " |r|=" << r << std::endl;
        }
    }
    double tcost = (clock() - start) / CLOCKS_PER_SEC;
    double omp_tcost = omp_get_wtime() - omp_start;
    std::cout << "itr=" << count << std::endl;
    std::cout << "Time cost = " << tcost << "(sec)\n";
    std::cout << "Time cost (WALL CLOCK)= " << omp_tcost << "(sec)\n";
    
    // check solution
    double r = 0.0;
    for (int i = 0; i < size; i++){
        double dr = sol[i] - xvec[i];
        r += dr * dr;
    }

    r = std::sqrt(r);
    std::cout << "|x - x0|=" << r << std::endl;
    return 0;
}
