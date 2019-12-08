//
// Large Scale Computing
// 2D Heat/Mass Transfer
// ex07.cpp : Numerical Analysis
// Solve for
// d2c/dx2 + d2c/dy2 + 1 = 0
// With the boundary conditions of c = 0 
// along lines of x=1,-1, and y=1, and -1
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

int main(int argc, char **argv){
  int num;
  std::cout << "Number of Points : ";
  std::cin >> num;

  if (std::cin.fail() || (num <= 3)) {
    std::cout << "BYE\n";
    return 0;
  }
  double *phi = new double[num * num];
  double *phi0 = new double[num * num];  

  double *rhs =  new double[num * num]; 
  double *ap =  new double[num * num]; 
  double *ae =  new double[num * num]; 
  double *aw =  new double[num * num];   
  double *an =  new double[num * num]; 
  double *as =  new double[num * num]; 

  int *ie =  new int[num * num]; 
  int *iw =  new int[num * num]; 
  int *in =  new int[num * num]; 
  int *is =  new int[num * num]; 

  /*assuming dx = dy : domain is 2x2 */
  double dx = 2.0 / (num - 1);

  /* Initialize */
  for (int n = 0 ; n < num * num ; n++){
    phi[n] = 0.0;
    phi0[n] = 0.0;
  }

  for (int n = 0 ; n < num * num ; n++){
    rhs[n] = -dx * dx;
    ap[n] = -4;
    ae[n] = 1;
    aw[n] = 1;
    an[n] = 1;
    as[n] = 1;
    ie[n] = n + 1;
    iw[n] = n - 1;
    in[n] = n + num;
    is[n] = n - num;
    if (((n % num) == 0) ||
	((n % num) == (num - 1)) ||
	((n / num) == 0) ||
	((n / num) == (num - 1))){
      rhs[n] = phi[n] * ap[n];
      ae[n] = aw[n] = an[n] = as[n] = 0;
      ie[n] = iw[n] = in[n] = is[n] = n;
    }
  }
  
  /* computing for phi with Jacobi Method */
  int k = 0;
  const double tol = 1.0e-8;
  double start = clock();
  while(1){
    for (int n = 0 ; n < num * num ; n++){
      phi[n] =  (rhs[n] - ae[n] * phi0[ie[n]] - aw[n] * phi0[iw[n]] - an[n] * phi0[in[n]] - as[n] * phi0[is[n]]) / ap[n];
    }
    double err = 0.0;
    for (int n = 0 ; n < num * num ; n++){
      double r = rhs[n] - ae[n] * phi[ie[n]] - aw[n] * phi[iw[n]] - an[n] * phi[in[n]] - as[n] * phi[is[n]] - ap[n] * phi[n];
      err += r * r;
      phi0[n] = phi[n];	
    }
    k++;
    if (err < tol * tol) break;
    if (k % 1000 == 0){
      std::cout << "# of Iteration=" << k << " err=" << std::sqrt(err) << std::endl;
    }
  }
  double tcost = (clock() - start) / CLOCKS_PER_SEC;
  std::cout << "# of Iteration=" << k << std::endl;
  std::cout << "Time cost = " << tcost << "(sec)\n";

  // Output Result
  std::ofstream ofile;
  ofile.open("ex07.dat");
  ofile << std::setprecision(16);
  ofile << std::scientific;
  for (int i = 0 ;i < num ; i++){
    for (int j = 0 ; j < num ; j++){
      double x = -1.0 + 2.0 * i / (num - 1);
      double y = -1.0 + 2.0 * j / (num - 1);
      ofile << x << " " << y << " " << phi[i + num * j] << std::endl;
    }
  }
  ofile.close();
  
  std::cout << "Done\n";
  return 0;
  printf("Done\n");
}
