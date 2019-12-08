//
//  Large Scale Computing
//  MPI Sample ex50.cpp
//
#include <iostream>
#include <cmath>
#include <iomanip>
#include <mpi.h>

int main(int argc, char *argv[]) {
	// Initialize
	MPI_Init(&argc,&argv);

	// get myid and # of processors 
	int numproc;
	int myid;
	MPI_Comm_size(MPI_COMM_WORLD,&numproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	if (myid == 0) {
		std::cout << "Isend/Irecev Test\n";
	}
	/* set values */
	int size = 10000;
	double *a = new double[size];
	double *b = new double[size];
	for (int i = 0; i < size; i++) {
		if (myid == 0) {
			a[i] = i;
			b[i] = size - i;
		} else {
			a[i] = 0;
			b[i] = 0;
		}
	}

	/* Send/Receive */
	MPI_Request r1;
	MPI_Request r2;
	if (myid == 0) {
		MPI_Isend(a, size, MPI_DOUBLE, 1, 21,MPI_COMM_WORLD,&r1);
		MPI_Isend(b, size, MPI_DOUBLE, 1, 22,MPI_COMM_WORLD,&r2);
	}
	if (myid == 1) {
		MPI_Irecv(b, size, MPI_DOUBLE, 0, 22,MPI_COMM_WORLD,&r2);
		MPI_Irecv(a, size, MPI_DOUBLE, 0, 21,MPI_COMM_WORLD,&r1);
	}

	// Can work to do while communicating. 

	// check if communication finished
	MPI_Status status;
	if (myid == 0) {
		MPI_Wait(&r1,&status);
		MPI_Wait(&r2,&status);
	}
	if (myid == 1) {
		MPI_Wait(&r2,&status);
		MPI_Wait(&r1,&status);
	}

	// wait until all processors come here 
	MPI_Barrier(MPI_COMM_WORLD);
	if (myid == 0) {
		std::cout << "Done\n";
	}
	MPI_Finalize();
	return 0;
}
