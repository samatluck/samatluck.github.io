//
//  Large Scale Computing
//  OpenMP Sample ex10.cpp
//

#include <iostream>
#include <omp.h>

int main (int argc, char *argv[]) {
  std::cout << "Start\n";

#pragma omp parallel
  {
    int id = omp_get_thread_num();
#pragma omp for
    for (int i = 0 ; i < 10 ; i++){
#pragma omp critical
      std::cout << "Thread" << id << " takes i=" << i << "\n";
    }
  }
  std::cout << "End\n";

  return 0;
}
