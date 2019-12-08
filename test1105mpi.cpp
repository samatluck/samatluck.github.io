//
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

int main(int argc, char **argv) {
	// Initialize
	MPI_Init(&argc, &argv);

	// get myid and # of processors
	int numproc;
	int myid;
	MPI_Comm_size(MPI_COMM_WORLD,&numproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	//Let num = 10000
	const int num = 10000;

	//Allocate two arrays of size num, LaTeX: A_i A i  and LaTeX: B_i B i .
	double *aArray = new double[num];
	double *bArray = new double[num];

	//Rank 0 assigns random values to those arrays.
	if (myid == 0) {
		for (int i = 0; i < num; i++) {
			aArray[i] = static_cast<double>(rand()) / RAND_MAX;
			bArray[i] = static_cast<double>(rand()) / RAND_MAX;
		}
	}

	//Rank 0 sends array elements to all other ranks.
	/* Broadcast */
	MPI_Bcast(aArray, num, MPI_DOUBLE, 0,MPI_COMM_WORLD);
	MPI_Bcast(bArray, num, MPI_DOUBLE, 0,MPI_COMM_WORLD);

	//Compute LaTeX: C_i=\sum_j^{num}{A_iB_j} C i = ∑ j n u m A i B j  in parallel.
	/* divide loop */
	int mystart = (num / numproc) * myid;
	int myend;
	if (num % numproc > myid) {
		mystart += myid;
		myend = mystart + (num / numproc) + 1;
	} else {
		mystart += num % numproc;
		myend = mystart + (num / numproc);
	}
	std::cout << "CPU" << myid << ":" << mystart << "~" << myend << std::endl;
	int mysize = myend - mystart;

	double *cArray = new double[mysize];

    for (int i = mystart; i < myend; i++) {
		cArray[i - mystart] = 0.0;
		for (int j = 0; j < num; j++) {
			cArray[i - mystart] += aArray[i] * bArray[j];
		}
	}

	//Compute |C| . Rank 0 gets the results.
	double cNorm = 0.0;
	for (int i = 0; i < mysize; i++) {
		cNorm += cArray[i] * cArray[i];
	}
	double tcNorm;
	MPI_Reduce(&cNorm, &tcNorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		tcNorm = std::sqrt(tcNorm);
		std::cout << "|C|=" << tcNorm << std::endl;
	}

	delete[] cArray;
	delete[] bArray;
	delete[] aArray;
	MPI_Finalize();
	return 0;
}
