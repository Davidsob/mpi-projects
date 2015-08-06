#include <mpi.h>

#ifndef PING_PONG_MAX
#define PING_PONG_MAX 10
#endif


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
    
    // play ping pong
    int count = 0;
    int partner_rank = (world_rank +1)%2; // 2 player!
    while (count < PING_PONG_MAX) {
        if(world_rank == count%2)
        {
            count++;
            MPI_Send(&count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
            printf("\nProcessor: %d sent and incremented ping-pong count "
                   "%d to partner: %d\n", world_rank, count, partner_rank);
        }else
        {
            count++;
            MPI_Recv(&count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("\nProcessor: %d received ping-pong count %d from partner: %d\n", world_rank, count, partner_rank);
        }
    }

    
    // Finalize the MPI environment
    MPI_Finalize();
}
