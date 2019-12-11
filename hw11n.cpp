/*
 * hw11.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: fuji
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <omp.h>

// Dimension
#define DIM 2
// Blob Size
#define EPSILON 0.005

// function declaration
double tsecond(); // timing method
double term1(double r, double ep);// function term 1
double term2(double r, double ep);// function term 2

// function definition
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
    double sq;
    sq = sqrt(r * r + ep * ep);
    return log(sq + ep) - ep * (sq + 2.0 * ep) / (sq + ep) / sq;
}

// function term 2
double term2(double r, double ep) {
    double sq;
    sq = sqrt(r * r + ep * ep);
    return (sq + 2.0 * ep) / (sq + ep) / (sq + ep) / sq;
}

// Main Routine
int main(int argc, char **argv) {
    

    
    double st = tsecond();
    const int numOfParticles = 1000;
    // Allocate space for position array
    double *loc = new double[numOfParticles * DIM];
    
    // Allocate space for force vector array
    double *foc = new double[numOfParticles * DIM];
    
    // Allocate space for velocity vector array
    double *vel = new double[numOfParticles * DIM];
    
    // Make Distribute particles and set forces
    for (int i = 0; i < numOfParticles; i++) {
        loc[i * DIM] = (double)rand() / RAND_MAX;
        loc[i * DIM + 1] = (double)rand() / RAND_MAX;
        foc[i * DIM] = (double)rand() / RAND_MAX - 0.5;
        foc[i * DIM + 1] = (double)rand() / RAND_MAX - 0.5;
    }
    
    const int num_dev = omp_get_num_devices();
    const int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    int totalProcs = omp_get_num_procs();
    
    /// Compute Velocities

#pragma omp target device(0)
#pragma omp parallel for
    for (int p = 0; p < numOfParticles; p++) {
        /* zeros */
        vel[p * DIM] = 0.0;
        vel[p * DIM + 1] = 0.0;
        
        /* loop for particles  */
#pragma omp parallel for
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

    vx /= numOfParticles;
    vy /= numOfParticles;
    
    // Show Results
    double et = tsecond();
    std::cout << "Mean Velocity = (" << vx << "," << vy << ")\n";
    std::cout << "Time cost = " << et - st << "(sec)\n";
    
    // cleanup
    delete [] loc;
    delete [] vel;
    delete [] foc;
    return 0;
}

