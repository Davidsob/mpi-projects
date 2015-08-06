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
    
    // hot potato
    int token = 0;
    if(world_rank != 0)
    {
        MPI_Recv(&token, 1, MPI_INT, world_rank-1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received token %d from process %d\n",world_rank,token,world_rank-1);
    }else{
        token = 13;
    }
    
    // now pass the message to your neighbor
    MPI_Send(&token, 1, MPI_INT,(world_rank+1)%world_size, 0, MPI_COMM_WORLD);
    
    if(world_rank == 0)
    {
        MPI_Recv(&token, 1, MPI_INT, world_size-1,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received token %d from process %d\n",world_rank,token,world_size-1);
    }

    
    // Finalize the MPI environment
    MPI_Finalize();
}
