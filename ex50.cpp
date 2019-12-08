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
	MPI::Init(argc, argv);

	// get myid and # of processors 
	int numproc = MPI::COMM_WORLD.Get_size();
	int myid = MPI::COMM_WORLD.Get_rank();
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
	MPI::Request r1;
	MPI::Request r2;
	if (myid == 0) {
		r1 = MPI::COMM_WORLD.Isend(a, size, MPI::INT, 1, 21);
		r2 = MPI::COMM_WORLD.Isend(b, size, MPI::INT, 1, 22);
	}
	if (myid == 1) {
		r2 = MPI::COMM_WORLD.Irecv(b, size, MPI::INT, 0, 22);
		r1 = MPI::COMM_WORLD.Irecv(a, size, MPI::INT, 0, 21);
	}

	// Can work to do while communicating. 

	// check if communication finished
	if (myid == 0) {
		r1.Wait();
		r2.Wait();
	}
	if (myid == 1) {
		r2.Wait();
		r1.Wait();
	}

	// wait until all processors come here 
	MPI::COMM_WORLD.Barrier();
	if (myid == 0) {
		std::cout << "Done\n";
	}
	MPI::Finalize();
	return 0;
}
