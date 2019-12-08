//
//  Large Scale Computing
//  MPI Sample ex49.cpp
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
	MPI::Status status;
	if (myid == 0) {
		MPI::COMM_WORLD.Send(a, size, MPI::DOUBLE, 1, 21);
		MPI::COMM_WORLD.Send(b, size, MPI::DOUBLE, 1, 22);
	}
	if (myid == 1) {
		MPI::COMM_WORLD.Recv(b, size, MPI::DOUBLE, 0, 22, status);
		MPI::COMM_WORLD.Recv(a, size, MPI::DOUBLE, 0, 21, status);
	}

	// wait until all processors come here 
	MPI::COMM_WORLD.Barrier();
	if (myid == 0) {
		std::cout << "Luckily came here\n";
	}
	MPI::Finalize();
	return 0;
}
