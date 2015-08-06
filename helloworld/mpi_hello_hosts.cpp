#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <mpi.h>

int main(int argc, char ** argv)
{
    MPI_Init(NULL,NULL);
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    char processor_name[512];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("Hello World from %s, proc %d out of %d total processors\n",processor_name,rank,nproc);

    MPI_Finalize();

}
