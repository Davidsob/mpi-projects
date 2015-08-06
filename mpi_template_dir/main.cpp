// This is a template for an mpi main.cpp file

#include <iostream>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char ** argv)
{
    MPI_Init(NULL,NULL);
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // insert code here
    fprintf(stdout,"Hello World: processor %d out of %d\n",rank,nproc);

    MPI_Finalize();
}
