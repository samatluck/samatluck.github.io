//
//  Large Scale Computing
//  OpenMP Sample ex11.cpp
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <omp.h>

int main (int argc, char *argv[]) {
  std::cout << "Start\n";

  double start = clock();
  double omp_start = omp_get_wtime();

#pragma omp parallel
  {
    int id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    if (id == 0){
      std::cout << nthreads << " threads work togather.\n";
    }
    // open file
    std::ostringstream fileName;
    fileName << "resEx11_" << std::setfill('0') << std::setw(2) << id << ".dat";
    std::ofstream ofile;
    ofile.open((fileName.str()).c_str());
    
    
    //
    //
    //
    //
    //#pragma omp for schedule(static,200)
    //#pragma omp for schedule(dynamic,500)
    //#pragma omp for schedule(guided,200)
    //#pragma omp for schedule(runtime)
#pragma omp for schedule(static,500)
    for (int i = 0 ; i < 10000 ; i++){
      for (int j = i ; j < 10000 ; j++){
	for (int k = 0 ; k < 100 ; k++){
	  /// weight
	}
      }
      ofile << i << " " << id << std::endl;
    }
    ofile.close();
  } // end of parallel
  double tcost = (clock() - start) / CLOCKS_PER_SEC;
  double omp_tcost = omp_get_wtime() - omp_start;
  std::cout << "Time cost (CPU) = " << tcost << "(sec)\n";
  std::cout << "Time cost (WALL CLOCK)= " << omp_tcost << "(sec)\n";
  std::cout << "End\n";

  return 0;
}
