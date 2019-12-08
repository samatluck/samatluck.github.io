// Eigen Test
// 
// ex30.cpp
// Need 'Eigen' library
#include <iostream>
#include <Eigen/Core>

int main(int argc, char **argv) {
  Eigen::Matrix2d mat1;
  mat1(0, 0) = 1.0;
  mat1(0, 1) = 2.0;
  mat1(1, 0) = 3.0;
  mat1(1, 1) = 4.0;

  Eigen::Matrix2d mat2;
  mat2(0, 0) = 3.0;
  mat2(0, 1) = 4.0;
  mat2(1, 0) = 5.0;
  mat2(1, 1) = 7.0;
	
  Eigen::Vector2d vec;
  vec(0) = 2.0;
  vec(1) = -3.0;

  // addition
  std::cout << "mat1 + mat2 = \n" << mat1 + mat2 << std::endl;
  // subtruction
  std::cout << "mat1 - mat2 = \n" << mat1 - mat2 << std::endl;
  // product
  std::cout << "mat1 * mat2 = \n" << mat1 * mat2 << std::endl;
  std::cout << "mat2 * mat1 = \n" << mat2 *  mat1 << std::endl;
  std::cout << "mat1 * vec  = \n" << mat1 * vec << std::endl;
  std::cout << "vec  * mat1 = \n" << vec.transpose() * mat1 << std::endl;
  // scalar multiplication
  std::cout << "mat1 *  5   = \n" << mat1 * 5.0 << std::endl;
  // scalar division
  std::cout << "mat1 /  5   = \n" << mat1 / 5.0 << std::endl;

  std::cout << "Here is the matrix mat1^T = \n" << mat1.transpose() <<std:: endl;

  return 0;
}
