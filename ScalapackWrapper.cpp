/*
 * ScalapackWrapper.cpp
 *
 *  Created on: Nov 13, 2013
 *      Author: fuji
 */
#include <mpi.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>

#include "ScalapackWrapper.h"

#ifdef MKL_ILP64
extern "C" {
  int numroc_(MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*, MKL_INT*);
  void descinit_(MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *);
  void pdgemm_(char *TRANSA, char *TRANSB, MKL_INT * M, MKL_INT * N, MKL_INT * K, double * ALPHA, double * A, MKL_INT * IA, MKL_INT * JA, MKL_INT * DESCA, double * B, MKL_INT * IB, MKL_INT * JB,
	       MKL_INT * DESCB, double * BETA, double * C, MKL_INT * IC, MKL_INT * JC, MKL_INT * DESCC);
}
#else

extern "C" {
  /* Cblacs declarations */
  void Cblacs_pinfo(int*, int*);
  void Cblacs_get(int, int, int*);
  void Cblacs_gridinit(int*, const char*, int, int);
  void Cblacs_pcoord(int, int, int*, int*);
  void Cblacs_gridexit(int);
  void Cblacs_barrier(int, const char*);
  void Cdgerv2d(int, int, int, double*, int, int, int);
  void Cdgesd2d(int, int, int, double*, int, int, int);

  int numroc_(int *, int*, int*, int*, int*);
  void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

  double pdlamch_(int *ictxt, char *cmach);
  double pdlange_(char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);

  void pdlacpy_(char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb);
  void pdgesv_(int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv, double *B, int *ib, int *jb, int *descb, int *info);
  void pdgemm_(char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA, double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB,
	       int * DESCB, double * BETA, double * C, int * IC, int * JC, int * DESCC);
  int indxg2p_(int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
}
#endif

ScalapackWrapper::ScalapackWrapper(int sizeOfMatrix, int sizeOfBlock, int sizeOfRowProcesses, int sizeOfColProcesses) :
  sizeOfMatrix(sizeOfMatrix), sizeOfBlock(sizeOfBlock), sizeOfRowProcesses(sizeOfRowProcesses), sizeOfColProcesses(sizeOfColProcesses) {
  /// All processes must call

  /// Initialize Process Grid
#ifdef MKL_ILP64
  MKL_INT nrhs = 1;
  MKL_INT izero = 0;
  MKL_INT ione = 1;
  MKL_INT nprow = sizeOfRowProcesses;
  MKL_INT npcol = sizeOfColProcesses;
  MKL_INT n = sizeOfMatrix;
  MKL_INT nb = sizeOfBlock;
  blacs_pinfo_(&myid, &numproc);
#else
  int nrhs = 1;
  int izero = 0;
  int ione = 1;
  int nprow = sizeOfRowProcesses;
  int npcol = sizeOfColProcesses;
  int n = sizeOfMatrix;
  int nb = sizeOfBlock;
  Cblacs_pinfo(&myid, &numproc);
#endif
  if (numproc != sizeOfRowProcesses * sizeOfColProcesses) {
    if (myid == 0) {
      std::cout << "use # of processors = " << sizeOfRowProcesses << "x" << sizeOfColProcesses << std::endl;
      std::cout << "You have " << numproc << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
#ifdef MKL_ILP64
  blacs_get_(&izero, &izero, &ictxt);
  blacs_gridinit_(&ictxt, "Row-major", &nprow, &npcol);
  blacs_pcoord_(&ictxt, &myid, &myrow, &mycol);
#else
  Cblacs_get(0, 0, &ictxt);
  Cblacs_gridinit(&ictxt, "Row-major", nprow, npcol);
  Cblacs_pcoord(ictxt, myid, &myrow, &mycol);
#endif
#ifdef DEBUG
  std::cout << "myid=" << myid << " mygrid=(" << myrow << "," << mycol << ")\n";
  std::cout << "n = " << sizeOfMatrix << "\tnrhs = " << nrhs << "\tprocess grid (" << nprow << "," << npcol << ")\t with blocks: (" << sizeOfBlock << "," << sizeOfBlock << ")\n";
#endif

  ///Compute the size of the local matrices (thanks to numroc)
  np = numroc_(&n, &nb, &myrow, &izero, &nprow);
  nq = numroc_(&n, &nb, &mycol, &izero, &npcol);
  nqrhs = numroc_(&nrhs, &nb, &mycol, &izero, &npcol);
#ifdef DEBUG
  std::cout << "myid=" << myid << " Local Matrix Size=(" << np << "," << nq << ")\n";
#endif
  /// Allocate and fill the matrices A and B
  matA = new double[np * nq];
  vecB = new double[np * nqrhs];
#ifdef MKL_ILP64
  ippiv = new MKL_INT[np + sizeOfBlock];
  MKL_INT itemp = np;
  MKL_INT info;
#else
  ippiv = new int[np + sizeOfBlock];
  int itemp = std::max(1, np);
  int info;
#endif
  ///    Initialize the array descriptor for the matrix A and B

  descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
  descinit_(descB, &n, &nrhs, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
}

ScalapackWrapper::~ScalapackWrapper() {
  /// All processes must call
  delete[] matA;
  delete[] vecB;
  delete[] ippiv;
}

double ScalapackWrapper::getValue(int gi, int gj) {
  /// local function, if (gi,gj) is not belong to myid, return 0;
  if (isMyOwn(gi, gj)) {
    int li = sizeOfBlock * (gi / (sizeOfBlock * sizeOfRowProcesses)) + gi % sizeOfBlock;
    int lj = sizeOfBlock * (gj / (sizeOfBlock * sizeOfColProcesses)) + gj % sizeOfBlock;
    return matA[lj + li * nq];
  }
  return 0;
}

bool ScalapackWrapper::setValue(int gi, int gj, double elm) {
  /// local function, if gi is not belong to myid, return false;
  if (isMyOwn(gi, gj)) {
    int li = sizeOfBlock * (gi / (sizeOfBlock * sizeOfRowProcesses)) + gi % sizeOfBlock;
    int lj = sizeOfBlock * (gj / (sizeOfBlock * sizeOfColProcesses)) + gj % sizeOfBlock;
    if ((li < 0) || (li >= np) || (lj < 0) || (lj >= nq)) {
      printf("out of range %d G(%d,%d) L(%d,%d)\n", myid, gi, gj, li, lj);
    }
    matA[lj + li * nq] = elm;
    return true;
  }
  return false;
}

double ScalapackWrapper::getValueRHS(int gi) {
  /// local function, if gi is not belong to myid, return 0;
  if (myrow == (gi / sizeOfBlock) % sizeOfRowProcesses) {
    int li = sizeOfBlock * (gi / (sizeOfBlock * sizeOfRowProcesses)) + gi % sizeOfBlock;
    return vecB[li * nqrhs];
  }
  return 0;
}

bool ScalapackWrapper::setValueRHS(int gi, double elm) {
  /// local function, if (gi,gi) is not belong to myid, return false;
  if (myrow == (gi / sizeOfBlock) % sizeOfRowProcesses) {
    int li = sizeOfBlock * (gi / (sizeOfBlock * sizeOfRowProcesses)) + gi % sizeOfBlock;
    vecB[li * nqrhs] = elm;
    return true;
  }
  return false;
}

bool ScalapackWrapper::isMyOwn(int gi, int gj) {
  /// local function, if (gi,gi) is not belong to myid, return false;
  int Pr = (gi / sizeOfBlock) % sizeOfRowProcesses;
  int Pc = (gj / sizeOfBlock) % sizeOfColProcesses;
  if ((Pr != myrow) || (Pc != mycol)) {
    return false;
  }
  return true;
}

void ScalapackWrapper::generateTestSystem() {
  int k = 0;
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < nq; j++) {
      int gi = sizeOfRowProcesses * sizeOfBlock * (i / sizeOfBlock) + myrow * sizeOfBlock + (i % sizeOfBlock);
      int gj = sizeOfColProcesses * sizeOfBlock * (j / sizeOfBlock) + mycol * sizeOfBlock + (j % sizeOfBlock);

      matA[k] = 1.0 / (double) (std::abs(gi - gj) + 10);
      k++;
    }
  }
  k = 0;
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < nqrhs; j++) {
      vecB[k] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
      k++;
    }
  }
}

void ScalapackWrapper::exportMat() {
  std::ostringstream fileName;
  fileName << "mat_" << std::setfill('0') << std::setw(2) << myid << ".dat";
  std::ofstream ofile;
  ofile.open((fileName.str()).c_str());
  int k = 0;
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < nq; j++) {
      ofile << matA[k] << " ";
      k++;
    }
    ofile << std::endl;
  }
  ofile << std::endl;
  ofile << std::endl;
  k = 0;
  for (int i = 0; i < np; i++) {
    for (int j = 0; j < nqrhs; j++) {
      ofile << matA[k] << " ";
      k++;
    }
    ofile << std::endl;
  }
  ofile << std::endl;
  ofile.close();
}

void ScalapackWrapper::pdgesv() {
  /// All processes must call
#ifdef MKL_ILP64
  MKL_INT nrhs = 1;
  MKL_INT izero = 0, ione = 1;
  MKL_INT info;
  MKL_INT nprow = sizeOfRowProcesses;
  MKL_INT npcol = sizeOfColProcesses;
  MKL_INT n = sizeOfMatrix;
#else
  int nrhs = 1;
  int izero = 0, ione = 1;
  int info;
  int nprow = sizeOfRowProcesses;
  int npcol = sizeOfColProcesses;
  int n = sizeOfMatrix;
#endif
#ifdef DEBUG
  ///Make a copy of A and the rhs for checking purposes
  double *matAcpy = new double[np * nq];
  double *vecX = new double[np * nqrhs];
  double *vecR = new double[np * nqrhs];
  pdlacpy_("All", &n, &n, matA, &ione, &ione, descA, matAcpy, &ione, &ione, descA);
  pdlacpy_("All", &n, &nrhs, vecB, &ione, &ione, descB, vecX, &ione, &ione, descB);
  if (myid == 0){
    printf("Solve with PDGESV \n");
  }
  double MPIt1 = MPI_Wtime();
#endif
  pdgesv_(&n, &nrhs, matA, &ione, &ione, descA, ippiv, vecB, &ione, &ione, descB, &info);

#ifdef DEBUG
  double MPIt2 = MPI_Wtime();
  std::cout << "INFO code returned by PDGESV = " << info << std::endl;
  ///		       Compute residual ||A * X  - B|| / ( ||X|| * ||A|| * eps * N )
  ///     Froebenius norm
  ///
  pdlacpy_("All", &n, &nrhs, vecX, &ione, &ione, descB, vecR, &ione, &ione, descB);
  double eps = pdlamch_(&ictxt, "Epsilon");
  double *work = new double;
  double mone=(-1.0e0),pone=(1.0e0);
  pdgemm_("N", "N", &n, &nrhs, &n, &pone, matAcpy, &ione, &ione, descA, vecB, &ione, &ione, descB, &mone, vecR, &ione, &ione, descB);
  double AnormF = pdlange_("F", &n, &n, matA, &ione, &ione, descA, work);
  double BnormF = pdlange_("F", &n, &nrhs, vecX, &ione, &ione, descB, work);
  double XnormF = pdlange_("F", &n, &nrhs, vecB, &ione, &ione, descB, work);
  double RnormF = pdlange_("F", &n, &nrhs, vecR, &ione, &ione, descB, work);
  double residF = RnormF / (AnormF * XnormF * eps * ((double) sizeOfMatrix));
  if (myid == 0) {
    std::cout << "\t||A * X  - B||_F / ( ||X||_F * ||A||_F * eps * N ) = " << residF << std::endl;
    if (residF < 10.0e+0){
      std::cout << "\tThe answer is correct.\n";
    }else{
      std::cout << "\tThe answer is suspicious.\n";
    }
    std::cout << "Time Cost for PDGESV = " << MPIt2 - MPIt1 << "(sec)\n";
  }
  delete [] matAcpy;
  delete [] vecX;
  delete [] vecR;
#endif

}
