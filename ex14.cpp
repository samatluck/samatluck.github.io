//
//  Large Scale Computing
//  OpenMP Sample ex14.cpp
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>

int main (int argc, char *argv[]) {
  int a = 0;
  double start = omp_get_wtime();
#pragma omp parallel
  {
    int id = omp_get_thread_num();

    // for loop without 'nowait'
#pragma omp single 
    {
      std::cout << "for loop\n";
    }

#pragma omp for
    for (int i = 0 ; i < 100; i++){
      for (int j = 0 ; j < 1000000 * (id + 1); j ++){
	// time weight
      }
    }
    double t1 = omp_get_wtime() - start;
#pragma omp critical (t1)
    {
      std::cout << "Thread" << id << " come here at " << t1 << std::endl;
    }

#pragma omp barrier
    // for loop with 'nowait'    
#pragma omp single 
    {
      std::cout << "for loop nowait\n";
    }
    start = omp_get_wtime();
#pragma omp for nowait
    for (int i = 0 ; i < 100; i++){
      for (int j = 0 ; j < 1000000 * (id + 1) ; j ++){
	// time weight
      }
    }
    double t2 = omp_get_wtime() - start;
#pragma omp critical (t2)
    {
      std::cout << "Thread" << id << " come here at " << t2 << std::endl;
    }

#pragma omp barrier
    // for loop with 'reduction' 'nowait'    
#pragma omp single 
    {
      std::cout << "for loop nowait reduction\n";
    }
    start = omp_get_wtime();
    a = 0;
#pragma omp for nowait reduction(+:a)
    for (int i = 0 ; i < 100 ; i++){
      a += 1; 
      for (int j = 0 ; j < 1000000 * (id + 1) ; j ++){
	// time weight
      }
    }
    double t3 = omp_get_wtime() - start;
#pragma omp critical (t3)
    {
      std::cout << "Thread" << id << " come here at " << t3 << " a= " << a << std::endl;
    }
#pragma omp barrier
#pragma omp critical (last)
    {
      std::cout << "Thread" << id << " a= " << a << std::endl;
    }
  }
  return 0;
}
