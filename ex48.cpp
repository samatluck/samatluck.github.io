//
//  Large Scale Computing
//  MPI Sample ex48.cpp
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

	/* set a value */
	int a;
	if (myid == 0) {
		a = 1234567;
	} else {
		a = 0;
	}

	std::cout << myid << " has a=" << a << std::endl;

	/* Send/Receive */
	MPI::Status status;
	if (myid == 0){
		MPI::COMM_WORLD.Send(&a, 1, MPI::INT, 1, 21);
	}
	if(myid == 1){
		MPI::COMM_WORLD.Recv(&a, 1, MPI::INT, 0, 21, status);
	}

	// wait until all processors come here 
	MPI::COMM_WORLD.Barrier();
	std::cout << "Now " << myid << " has a=" << a << std::endl;

	MPI::Finalize();
	return 0;
}
