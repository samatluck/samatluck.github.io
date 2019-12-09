//
//  Large Scale Computing
//  2D Heat/Mass Transfer
//  ex57.cpp
//  Solve for
//  d2c/dx2 + d2c/dy2 + 1 = 0
//  With the boundary conditions of c = 0
//  along lines of x=1,-1, and y=1, and -1
//  use PETSc
//

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <petsc.h>

#define NDIM 2

int main(int argc, char **argv) {
#if PETSC_VERSION_LE(3,1,0)
	PetscTruth fl1,fl2,fl3;
#else
	PetscBool fl1, fl2, fl3;
#endif
	/* Initialize PETSc and MPI */
	PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);

	/* set parameters */
	int p, q, num;
#if PETSC_VERSION_LE(3,6,0)	
	PetscOptionsGetInt(PETSC_NULL, "-n", &num, &fl1);
	PetscOptionsGetInt(PETSC_NULL, "-p", &p, &fl2);
	PetscOptionsGetInt(PETSC_NULL, "-q", &q, &fl3);
#else
	PetscOptionsGetInt(NULL, NULL, "-n", &num, &fl1);
	PetscOptionsGetInt(NULL, NULL, "-p", &p, &fl2);
	PetscOptionsGetInt(NULL, NULL, "-q", &q, &fl3);
#endif	
	if ((fl1 == PETSC_FALSE) || (fl2 == PETSC_FALSE) || (fl3 == PETSC_FALSE)) {
		PetscPrintf(PETSC_COMM_WORLD, "Usage:%s -n [NUM] -p [NProw] -q [NPcol] \n", argv[0]);
		PetscFinalize();
		return 0;
	}

	/*assuming dx = dy : domain is 2x2 */
	double dx = 2.0 / (double) (num - 1);

	/* get # of process and myid, use MPI commands */
	PetscMPIInt numproc, myid;
	MPI_Comm_size(PETSC_COMM_WORLD, &numproc);
	MPI_Comm_rank(PETSC_COMM_WORLD, &myid);

	if (p * q != numproc) {
		PetscPrintf(PETSC_COMM_WORLD, "[NProw] * [NPcol] must be %d\n", numproc);
		PetscFinalize();
		return 0;
	}

	/* Create Processor Grid */
	int dim[NDIM], period[NDIM], reorder;
	int coord[NDIM];

	dim[0] = p; /* p x q grid */
	dim[1] = q;
	period[0] = 0; /* no periodic */
	period[1] = 0; /* no periodic */
	reorder = 1; /* reorder = yes */
	MPI_Comm new_comm;
	MPI_Cart_create(PETSC_COMM_WORLD, NDIM, dim, period, reorder, &new_comm);
	MPI_Cart_coords(new_comm, myid, NDIM, coord);

	/* start end */
	int myXstart = (num / p) * coord[0] + ((num % p) < coord[0] ? (num % p) : coord[0]);
	int myXend = myXstart + (num / p) + ((num % p) > coord[0]);
	int myXsize = myXend - myXstart;

	int myYstart = (num / q) * coord[1] + ((num % q) < coord[1] ? (num % q) : coord[1]);
	int myYend = myYstart + (num / q) + ((num % q) > coord[1]);
	int myYsize = myYend - myYstart;

	int mysize = myXsize * myYsize;
	PetscSynchronizedPrintf(new_comm, "Proc%d: X %d~%d Y %d~%d\n", myid, myXstart, myXend, myYstart, myYend);
#if PETSC_VERSION_LE(3,5,0)
  PetscSynchronizedFlush(new_comm);
#else
	PetscSynchronizedFlush(new_comm, PETSC_STDOUT);
#endif

	/* define vector */
	Vec phi;
	Vec rhs;
	VecCreateMPI(new_comm, mysize, num * num, &phi);
	VecDuplicate(phi, &rhs);
	PetscInt mystart, myend;
	VecGetOwnershipRange(rhs, &mystart, &myend);
	PetscSynchronizedPrintf(new_comm, "Proc%d: Vec phi %d~%d\n", myid, mystart, myend);
#if PETSC_VERSION_LE(3,5,0)
  PetscSynchronizedFlush(new_comm);
#else
	PetscSynchronizedFlush(new_comm, PETSC_STDOUT);
#endif

	/* define matrix */
	Mat A;
#if PETSC_VERSION_LE(3,2,0)
	MatCreateMPIAIJ(new_comm,
#else
	MatCreateAIJ(new_comm,
#endif
			myend - mystart, myend - mystart, // Local Size
			num * num, num * num, // Global Size
			PETSC_DECIDE, NULL, /* number of non-zero (PETSc to decide)*/
			PETSC_DECIDE, NULL, &A);

	/* make domain grid indices -> vector index */
	int *gIndexBuf = new int[num * num];
	for (int i = 0; i < num * num; i++) {
		gIndexBuf[i] = 0;
	}
	int **gIndex = new int*[num];
	for (int i = 0; i < num; i++)
		gIndex[i] = &gIndexBuf[i * num];
	int n = mystart;
	for (int j = myYstart; j < myYend; j++) {
		for (int i = myXstart; i < myXend; i++) {
			gIndex[i][j] = n;
			n++;
		}
	}
	int *gIndexBuf2 = new int[num * num];
	MPI_Allreduce(gIndexBuf, gIndexBuf2, num * num, MPI_INT, MPI_SUM, new_comm);
	for (int i = 0; i < num * num; i++) {
		gIndexBuf[i] = gIndexBuf2[i];
	}
	delete[] gIndexBuf2;

	/* set matrix & rhs */
	for (int j = myYstart; j < myYend; j++) {
		for (int i = myXstart; i < myXend; i++) {

			int ip = gIndex[i][j];
			double ae = 1;
			double aw = 1;
			double an = 1;
			double as = 1;
			double ap = ae + aw + an + as;

			if (i == 0) { /* left */
				MatSetValue(A, ip, ip, 1.0, INSERT_VALUES);
				VecSetValue(rhs, ip, 0.0, INSERT_VALUES);
				continue;
			}
			if (i == (num - 1)) { /* right */
				MatSetValue(A, ip, ip, 1.0, INSERT_VALUES);
				VecSetValue(rhs, ip, 0.0, INSERT_VALUES);
				continue;
			}
			if (j == 0) { /* bottom */
				MatSetValue(A, ip, ip, 1.0, INSERT_VALUES);
				VecSetValue(rhs, ip, 0.0, INSERT_VALUES);
				continue;
			}
			if (j == (num - 1)) { /* top */
				MatSetValue(A, ip, ip, 1.0, INSERT_VALUES);
				VecSetValue(rhs, ip, 0.0, INSERT_VALUES);
				continue;
			}

			int ie = gIndex[i + 1][j];
			int iw = gIndex[i - 1][j];
			int in = gIndex[i][j + 1];
			int is = gIndex[i][j - 1];
			MatSetValue(A, ip, ip, -ap, INSERT_VALUES);
			MatSetValue(A, ip, ie, ae, INSERT_VALUES);
			MatSetValue(A, ip, iw, aw, INSERT_VALUES);
			MatSetValue(A, ip, in, an, INSERT_VALUES);
			MatSetValue(A, ip, is, as, INSERT_VALUES);
			VecSetValue(rhs, ip, -dx * dx, INSERT_VALUES);
		}
	}

	/* assemble */
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

	/* SOLVE LINEAR SYSTEM */
	PetscLogDouble kspBegin, kspEnd;
	kspBegin = MPI_Wtime();
	KSP ksp; /* linear solver context */
	PC pc; /* preconditioner context */
	KSPCreate(new_comm, &ksp); /* create KSP object */
#if PETSC_VERSION_LE(3,5,0)
	KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);
#else
	KSPSetOperators(ksp, A, A);
#endif
	KSPGetPC(ksp, &pc); /* create pre-conditionar object */
	KSPSetFromOptions(ksp);
	KSPSolve(ksp, rhs, phi);
	kspEnd = MPI_Wtime();
	KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD);
	PetscPrintf(new_comm, "Time Cost = %e(s)\n", kspEnd - kspBegin);

	/* Output Result */
	PetscScalar *pphi;
	std::ostringstream fileName;
	fileName << "resEx57_" << std::setfill('0') << std::setw(2) << myid << ".dat";
	std::ofstream ofile;
	ofile.open((fileName.str()).c_str());
	ofile << std::setprecision(16);
	ofile << std::scientific;
	VecGetArray(phi, &pphi);
	n = 0;
	for (int j = myYstart; j < myYend; j++) {
		for (int i = myXstart; i < myXend; i++) {
			double x = -1.0 + 2.0 * (double) i / (double) (num - 1);
			double y = -1.0 + 2.0 * (double) j / (double) (num - 1);
			ofile << x << " " << y << " " << pphi[n] << std::endl;
			n++;
		}
	}
	VecRestoreArray(phi, &pphi);
	ofile.close();

#if PETSC_VERSION_LE(3,1,0)
	VecDestroy(phi);
	VecDestroy(rhs);
	MatDestroy(A);
	KSPDestroy(ksp);
#else
	VecDestroy(&phi);
	VecDestroy(&rhs);
	MatDestroy(&A);
	KSPDestroy(&ksp);
#endif
	PetscFinalize();
	return 0;
}
