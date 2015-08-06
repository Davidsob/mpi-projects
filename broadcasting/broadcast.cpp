#include <mpi.h>


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
  MPI_Init(NULL,NULL);
  
  int world_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&world_size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // use my broadcast
  int number;
  if(rank == 0)
  {
    number = 13;
    myBcast(&number,1,MPI_INT,0,MPI_COMM_WORLD);
  }else{
    myBcast(&number,1,MPI_INT,0,MPI_COMM_WORLD);
    printf("Process %d received data %d from process %d\n",rank,0,number);
  }
  
  MPI_Finalize();
}
