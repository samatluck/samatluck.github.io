/*
  Large Scale Computing
  1D Heat/Mass Transfer
  ex50.c : use PETSc
 */
#include "petscksp.h"
PetscLogDouble t_cost();
int main(int argc, char **argv){
  double a,b,d;
  Vec phi;
  Vec rhs;
  Mat A;
  KSP ksp;         /* linear solver context */
  PC pc;           /* preconditioner context */
  PetscScalar *pphi;
  double dx;
  PetscInt mystart,myend;
  PetscMPIInt numproc,myid;
  int i,proc,num;
  FILE *fp;

  /* Initilize PETSc and MPI */
  PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

  /* set parameters */
  PetscOptionsGetInt(PETSC_NULL,"-n",&num,PETSC_NULL);
  PetscOptionsGetScalar(PETSC_NULL,"-a",&a,PETSC_NULL);
  PetscOptionsGetScalar(PETSC_NULL,"-b",&b,PETSC_NULL);
  PetscOptionsGetScalar(PETSC_NULL,"-d",&d,PETSC_NULL);
  
  /* get # of process and myid, use MPI commands */
  MPI_Comm_size(PETSC_COMM_WORLD,&numproc);
  MPI_Comm_rank(PETSC_COMM_WORLD,&myid);

  if (argc < 5){
    PetscPrintf(PETSC_COMM_WORLD,"Usage:%s -n [NUM] -a [A] -b [B] -d [D]\n",argv[0]);
    MPI_Abort(PETSC_COMM_WORLD, -1);
  }
  
  /* set grid size */
  dx = 1.0 / (double)(num -1);
  printf("num=%d A=%e B=%e D=%e dx=%e\n",num,a,b,d,dx);

  /* define vector */
  VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,num,&phi);
  VecDuplicate(phi,&rhs); 
  VecGetOwnershipRange(rhs,&mystart,&myend);  

  /* define matrix */
#if (PETSC_VERSION_MAJOR >= 3) && (PETSC_VERSION_MINOR >= 3)
  MatCreateAIJ(PETSC_COMM_WORLD,
#else
  MatCreateMPIAIJ(PETSC_COMM_WORLD,
#endif
		  myend - mystart,myend - mystart,
		  num,num,      /* size of matrix */
		  0,PETSC_NULL, /* number of non-zero (no idea)*/
		  0,PETSC_NULL, /* number of non-zero (no idea)*/
		  &A);  
  
  /* set matrix rhs */
  for (i = mystart ; i < myend ; i++){
    if (i == 0){
      VecSetValue(rhs,0, a, INSERT_VALUES);    
      MatSetValue(A, 0, 0, 1.0, INSERT_VALUES);       
      continue;
    }
    if (i == (num - 1)){
      VecSetValue(rhs,num-1, b, INSERT_VALUES); 
      MatSetValue(A, num-1, num-1, 1.0, INSERT_VALUES);            
      continue;
    }
    VecSetValue(rhs,i, -dx * dx / d * (dx * (double)i), INSERT_VALUES);     
    MatSetValue(A, i, i - 1, 1.0, INSERT_VALUES);   
    MatSetValue(A, i, i,    -2.0, INSERT_VALUES);   
    MatSetValue(A, i, i + 1, 1.0, INSERT_VALUES);   
  }

  /* assemble */
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);

  t_cost();
  KSPCreate(PETSC_COMM_WORLD,&ksp); /* create KSP object */
  KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
  KSPGetPC(ksp,&pc);  /* create pre-conditionar object */
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,rhs,phi);
  PetscPrintf(PETSC_COMM_WORLD,"time = %e sec\n",t_cost());

  /* 
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
 
  /* Output Result */
  for (proc = 0 ; proc < numproc ; proc++){
    if (myid == proc){
      if (myid == 0) fp = fopen("res.dat","w");
      else fp = fopen("res.dat","a");

      VecGetArray(phi, &pphi);
      for (i = mystart ; i < myend ; i++){
	fprintf(fp,"%e %e\n",dx * (double)i,pphi[i - mystart]);
      }
      VecRestoreArray(phi, &pphi);
      fclose(fp);
    }
    PetscBarrier(PETSC_NULL);
  }
  
  VecDestroy(phi);
  VecDestroy(rhs);
  MatDestroy(A);
  KSPDestroy(ksp);

  PetscFinalize();
  return 0;
}

/* for timing */
PetscLogDouble t_cost(){
  static PetscLogDouble v1 = 0.0;
  PetscLogDouble v2,et;
  PetscGetTime(&v2);
  et = v2 - v1;
  v1 = v2;
  return(et);
} 
