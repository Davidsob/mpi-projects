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
    
    // get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len = 0;
    MPI_Get_processor_name(processor_name, &name_len);
    
    // print a hello world message
    printf("Hello world from processor %s, rank %d"
           "    out of %d processors.\n",processor_name,world_rank, world_size);
    
    // Finalize the MPI environment
    MPI_Finalize();
}
