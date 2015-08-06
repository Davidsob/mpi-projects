#include <mpi.h>

int main(int argc, char ** argv)
{
    // initialize MPI environment
    MPI_Init(NULL,NULL);
    
    
    // get the number of processes
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    
    // get rank of the process
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // play send and recieve between to processors
    int number = 0;
    if(world_rank == 0){
        number = -1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }else if(world_rank == 1)
    {
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1 recieved number %d from process 0\n",number);
    }
    
    // Finalize the MPI environment
    MPI_Finalize();
}
