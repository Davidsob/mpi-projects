#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
void myBcast(void * data, int count, MPI_Datatype datatype,int root,MPI_Comm comm)
{
  int world_size, rank;
  MPI_Comm_size(comm,&world_size);
  MPI_Comm_rank(comm,&rank);
  
  if(rank == root)
  {
    int i;
    for(i = 0; i < world_size; i++)
    {
      if(i != root)
      {
        MPI_Send(data,count,datatype,i,0,comm);
      }
    }
  }else{
    MPI_Recv(data,count,datatype,root,0,comm,MPI_STATUS_IGNORE);
  }
}

int main(int argc, char ** argv)
{
  if(argc != 3){
    fprintf(stderr,"Usage: compare_bcast num_elements num_trials\n");
    exit(1);
  }
  
  int num_elements = atoi(argv[1]);
  int num_trials = atoi(argv[2]);
  
  MPI_Init(NULL,NULL);
  
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  double my_bcast_time = 0.0;
  double bcast_time = 0.0;
  int * data = (int *)calloc(num_elements,sizeof(int));
  // use my broadcast
  for(int i = 0; i < num_trials; i++)
  {
    MPI_Barrier(MPI_COMM_WORLD); // sync before starting timer
    my_bcast_time -= MPI_Wtime();
    myBcast(data,num_elements,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // sync before getting final time
    my_bcast_time += MPI_Wtime();
    
    // time mpi tree based broadcast
    MPI_Barrier(MPI_COMM_WORLD); // sync before starting timer
    bcast_time -= MPI_Wtime();
    MPI_Bcast(data,num_elements,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); // sync before getting final time
    bcast_time += MPI_Wtime();
  }
  
  // print timing
  if(rank == 0)
  {
    printf("Data size = %d bytes, Trials = %d\n", num_elements * (int)sizeof(int),
           num_trials);
    printf("Avg my_bcast time = %lf\n", my_bcast_time / num_trials);
    printf("Avg MPI_Bcast time = %lf\n", bcast_time / num_trials);
  }
  
  free(data);
  MPI_Finalize();
}
