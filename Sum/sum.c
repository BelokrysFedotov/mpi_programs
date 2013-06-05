#include <stdio.h>
#include <mpi.h>

#define MIN 1
#define MAX 1000000
#define STEP 6

int main(int argc, char *argv[]) {
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  double curtime;

  long a,b,i,j,k,sum,ssum;



  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &namelen);

  printf("Process %d on %s out of %d\n", rank, processor_name, numprocs);

  a = MIN+(MAX-MIN)/numprocs * rank;
  b = MIN+(MAX-MIN)/numprocs * (rank+1) -1;
  if( rank == 0) a = MIN;
  if( rank+1 == numprocs ) b = MAX;
  
  i = MIN%STEP;  
  k = a%STEP;
  j = a;
  while(k!=i){j++;k++;k=k%STEP;}

  sum=0;
  for(i=j;i<=b;i+=STEP){
    sum+=i;
  }
 
  printf("p %d: a=%ld\tb=%ld\tsum=%ld\n",rank,a,b,sum);

  MPI_Reduce(&sum,&ssum,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);

  if(rank==0){
    printf("ssum=%ld\n",ssum);
  }

  MPI_Finalize();
}
