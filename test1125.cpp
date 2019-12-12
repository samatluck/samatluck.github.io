// Practice MIC offloading
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>

int main(int argc, char **argv){
    //Let num = 10000
    const int num = 10000;
    
    int num_dev = omp_get_num_devices();
    std::cout << "number of devices " << num_dev << std::endl;
    
    //Allocate two arrays of size num, A_i and B_i
    double *aArray = new double[num];
    double *bArray = new double[num];
    
    //Rank 0 assigns random values to those arrays.
    for (int i = 0 ; i < num ; i++){
        aArray[i] = static_cast<double>(rand()) / RAND_MAX;
        bArray[i] = static_cast<double>(rand()) / RAND_MAX;
    }
    
    //Allocate array of size num, C_i
    double *cArray = new double[num];
    
    int numproc = num_dev + 1;
#pragma omp parallel num_threads(numproc)
#pragma omp single
    {
        for (int dev = 0; dev < numproc ; dev++) {
#pragma omp task firstprivate(dev)
            {
                /* divide domain */
                int mystart = (num / numproc) * dev;
                int myend;
                if (num % numproc > dev) {
                    mystart += dev;
                    myend = mystart + (num / numproc) + 1;
                } else {
                    mystart += num % numproc;
                    myend = mystart + (num / numproc);
                }
                int mysize = myend - mystart;
                
                // allocate local cArray
                double *cArray_dev = new double[mysize];
#pragma omp target if(dev != num_dev) device(dev) map(to:aArray[mystart:myend]) map(to:bArray[0:num]) map(from:cArray_dev[0:mysize])
                {// offload begins Transfer aArray[mystart:myend] bArray[0:num] from host to device.
#pragma omp parallel for
                    for (int i = mystart ; i < myend ; i++){
                        cArray_dev[i - mystart] = 0.0;
                        for (int j = 0 ; j < num ; j++){
                            cArray_dev[i - mystart] += aArray[i] * bArray[j];
                        }
                    }
                }//  offload ends Transfer cArray_dev from device to host.
                for (int i = mystart ; i < myend ; i++){
                    cArray[i] = cArray_dev[i - mystart];
                }
                delete [] cArray_dev;
            }
        }
    }
    
    //Compute ||C|| . Host gets the results.
    double cNorm = 0.0;
    for (int i = 0 ; i < num ; i++){
        cNorm += cArray[i] * cArray[i];
    }
    cNorm = std::sqrt(cNorm);
    std::cout << "||C||=" << cNorm << std::endl;
    
    delete [] cArray;
    delete [] bArray;
    delete [] aArray;
    return 0;
}
