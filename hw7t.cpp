#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <omp.h>
#include <iomanip>
#include <iomanip>
#include <mpi.h>

// Dimension
#define DIM 2
// Blob Size
#define EPSILON 0.005

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

// function term 1
double term1(double r, double ep) {
#ifdef VALGRIND
    return 1.0;
#else
    double sq;
    
    sq = sqrt(r * r + ep * ep);
    
    return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
#endif
}

// function term 2
double term2(double r, double ep) {
#ifdef VALGRIND
    return 1.0;
#else
    double sq;
    
    sq = sqrt(r * r + ep * ep);
    
    return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
#endif
}

// Main Routine
int main(int argc, char **argv) {
    
    MPI_Init(&argc, &argv);
    int numproc;
    int myid;
    MPI_Comm_size(MPI_COMM_WORLD,&numproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    
    double st = tsecond();
    const int numOfParticles = 100000;
    // Allocate space for position array
    double *loc = new double[numOfParticles * DIM];
    
    // Allocate space for force vector array
    double *foc = new double[numOfParticles * DIM];
    
    // Allocate space for velocity vector array
    double *vel = new double[numOfParticles * DIM];
    
    // Make Distribute particles and set forces
    if (myid == 0) {
        for (int i = 0; i < numOfParticles; i++) {
            loc[i * DIM] = (double)rand() / RAND_MAX;
            loc[i * DIM + 1] = (double)rand() / RAND_MAX;
            foc[i * DIM] = (double)rand() / RAND_MAX - 0.5;
            foc[i * DIM + 1] = (double)rand() / RAND_MAX - 0.5;
        }
    }
    //Rank 0 sends array elements to all other ranks.
    /* Broadcast */
    MPI_Bcast(loc, numOfParticles * DIM, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(foc, numOfParticles * DIM, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    
  
   
    /// Compute Velocities
    for (int p = 0; p < numOfParticles; p++) {
        /* zeros */
        vel[p * DIM] = 0.0;
        vel[p * DIM + 1] = 0.0;
        
        /* loop for particles  */
        for (int i = 0; i < numOfParticles; i++) {
            double dx = loc[p * DIM] - loc[i * DIM];
            double dy = loc[p * DIM + 1] - loc[i * DIM + 1];
            double r = sqrt(dx * dx + dy * dy);
            
            double tr1 = term1(r, EPSILON) / (4.0 * M_PI);
            double tr2 = term2(r, EPSILON) / (4.0 * M_PI);
            
            tr2 *= foc[i * DIM] * dx + foc[i * DIM + 1] * dy;
            
            vel[p * DIM] += -foc[i * DIM] * tr1 + tr2 * dx;
            vel[p * DIM + 1] += -foc[i * DIM + 1] * tr1 + tr2 * dy;
        }
    }
    
    // Compute Average Velocity
    double vx = 0.0;
    double vy = 0.0;
    for (int i = 0; i < numOfParticles; i++) {
        vx += vel[i * DIM];
        vy += vel[i * DIM + 1];
    }
    
   // MPI_Reduce(&vx,&vy,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    vx /= numOfParticles;
    vy /= numOfParticles;
    
    // Show Results
    double et = tsecond();
    if (myid == 0) {
        std::cout << "Mean Velocity = (" << vx << "," << vy << ")\n";
    }
    std::cout << "Time cost = " << et - st << "(sec)\n";
    
    // cleanup
    delete [] loc;
    delete [] vel;
    delete [] foc;
    MPI_Finalize();
    return 0;
}

