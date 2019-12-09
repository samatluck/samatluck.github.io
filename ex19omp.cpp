//
// Large Scale Computing
// 3D Heat/Mass Transfer
// ex19c.cpp :
// Solve for
// d2c/dx2 + d2c/dy2 + d2c/dz2 + 1 = 0
// With the boundary conditions of c = 0
// along planes of x=1,-1,  y=1,-1, and z=1,-1
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <omp.h>

/// Allocate 3D array
double ***makeDoubleArray3D(double *&v,int s1,int s2,int s3){
  v = new double[s1 * s2 * s3];
  double ***a = new double**[s1];
  for(int i = 0; i < s1; ++i){
    a[i] = new double*[s2];
    for(int j = 0; j < s2; ++j){
	a[i][j] = v + s2 * s3 * i + s3 * j;
    }
  }
  return a;
}

/// Destroy 3D array
void destroyDoubleArray3D(double *v,double ***m,int s1,int s2,int s3){
  for(int i = 0; i < s1; ++i){
    delete [] m[i];
  }
  delete [] m;
  delete [] v;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << argv[0] << " [number of points]\n";
    return 0;
  }

  int num = std::atoi(argv[1]);
  std::cout << "Number of Points=" << num << std::endl;

  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0){
      std::cout << "Number of Threads = " << omp_get_num_threads() << std::endl; 
    }
  }
  double *phi;
  double *phi0;
  double ***phi3D = makeDoubleArray3D(phi,num,num,num);
  double ***phi03D = makeDoubleArray3D(phi0,num,num,num);


  /*assuming dx = dy = dz: domain is 2x2x2 */
  double dx = 2.0 / (num - 1);

  /* Initialize */
  std::fill(phi,phi+num * num * num,0.0);
  std::fill(phi0,phi0+num * num * num,0.0);

  /* computing for phi with Jacobi Method */
  int itc = 0;
  const double tol = 1.0e-8;
  double start = clock();
  double omp_start = omp_get_wtime();
  while (1) {

#pragma omp parallel for collapse(3)
    for (int i = 1 ; i < num - 1 ; i++){	  
      for (int j = 1 ; j < num - 1 ; j++){
	for (int k = 1 ; k < num - 1 ; k++){
	  double sum = dx * dx;
	  sum -= phi03D[i-1][j][k];
	  sum -= phi03D[i+1][j][k];
	  sum -= phi03D[i][j-1][k];
	  sum -= phi03D[i][j+1][k];
	  sum -= phi03D[i][j][k-1];
	  sum -= phi03D[i][j][k+1];		
	  phi3D[i][j][k] = sum / (-6.0);
	}
      }
    }

    // compute residual 
    double err = 0.0;
#pragma omp parallel for collapse(3) reduction(+:err)    
    for (int i = 1 ; i < num - 1 ; i++){	  
      for (int j = 1 ; j < num - 1 ; j++){
	for (int k = 1 ; k < num - 1 ; k++){
	  double r = dx * dx;
	  r -= phi3D[i-1][j][k];
	  r -= phi3D[i+1][j][k];
	  r -= phi3D[i][j-1][k];
	  r -= phi3D[i][j+1][k];
	  r -= phi3D[i][j][k-1];
	  r -= phi3D[i][j][k+1];		
	  r += 6 * phi3D[i][j][k];
	  err += r * r;
	  phi03D[i][j][k] = phi3D[i][j][k];
	}
      }
    }

    itc++;
    if (err < tol * tol)
      break;
    if (itc % 1000 == 0) {
      std::cout << "# of Iteration=" << itc
		<< " err=" << std::sqrt(err) << std::endl;
    }
  }
  double tcost = (clock() - start) / CLOCKS_PER_SEC;
  double omp_tcost = (omp_get_wtime() - omp_start);
  std::cout << "Number of Iteration=" << itc << std::endl;
  std::cout << "Time cost (CPU) = " << tcost << "(sec)\n";
  std::cout << "Time cost (WTIME) = " << omp_tcost << "(sec)\n";

  // Output Result
  std::ofstream ofile;
  ofile.open("ex19c.dat");
  ofile << std::setprecision(16);
  ofile << std::scientific;
  for (int i = 0; i < num; i++) {
    for (int j = 0; j < num; j++) {
      for (int k = 0; k < num; k++) {
	double x = -1.0 + 2.0 * i / (num - 1);
	double y = -1.0 + 2.0 * j / (num - 1);
	double z = -1.0 + 2.0 * k / (num - 1);
	ofile << x << " " << y << " " << z << " "
	      << phi3D[i][j][k] << std::endl;
      }
    }
  }
  ofile.close();

  // cleanup
  destroyDoubleArray3D(phi,phi3D,num,num,num);
  destroyDoubleArray3D(phi0,phi03D,num,num,num);

  std::cout << "Done\n";
  return 0;
}

