
//  Large Scale Computing
//  Homework
//
#include <iostream>
#include <cmath>
#include <cstdlib>
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
    
    /* initialize random seed: */
    srand (myid);
    
    /* set values */
    const int size = 10000;
    double *aArray = new double[size];
    double *bArray = new double[size];
    
    // compute average
    double average = 0;
    for (int i = 0; i < size; i++) {
        aArray[i] = static_cast<double>(rand()) / RAND_MAX;
        average += aArray[i];
        bArray[i] = 0;
    }
    average /= size;
    
    // setup to/from ranks and tags
    int rightProc = (myid + 1) % numproc;
    int leftProc = (myid - 1 + numproc) % numproc;
    int tagSend = myid + rightProc * numproc;
    int tagRecv = leftProc + myid * numproc;
    
    // print in a right order
    for (int p = 0 ; p < numproc ; p++){
        MPI_Barrier(MPI_COMM_WORLD);
        if (p == myid){
            std::cout << "Rank" << myid << " sends data of average=" << average << " to Rank" << rightProc << " with tag=" << tagSend << std::endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0){
        std::cout << "Data send/receive\n";
    }
    
    MPI::Request r1;
    MPI::Request r2;
   
    r1 = MPI::COMM_WORLD.Isend(aArray, size, MPI_DOUBLE, leftProc, tagSend); //,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    r2 = MPI::COMM_WORLD.Irecv(bArray, size, MPI_DOUBLE, rightProc, tagRecv);
    MPI_Barrier(MPI_COMM_WORLD);
   
    // compute average
    average = 0;
    for (int i = 0; i < size; i++) {
        average += bArray[i];
    }
    average /= size;
    
    // print in a right order
    for (int p = 0 ; p < numproc ; p++){
        MPI_Barrier(MPI_COMM_WORLD);
        if (p == myid){
            std::cout << "Rank" << myid << " received data of average=" << average << " from Rank" << leftProc << " with tag=" << tagRecv << std::endl;
        }
    }
    
    // wait until all processors come here
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid == 0) {
        std::cout << "Done\n";
    }
    MPI_Finalize();
    return 0;
}
