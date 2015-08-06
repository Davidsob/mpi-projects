#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <Random.h>

double *createRand(int num_els)
{
    Random rg;
    double * nums = new double[num_els];
    for(int i = 0; i < num_els; i++)
      nums[i] = rg.getNormal(13,4);

    return nums;
}


double computeAvg(double * nums,int num_els)
{
    double sum = 0.0;
    for(int i = 0; i < num_els; i++)
      sum += nums[i];

    return sum/num_els;
}

int main(int argc, char ** argv)
{
    MPI_Init(NULL,NULL);
    int rank, world;
    MPI_Comm_size(MPI_COMM_WORLD,&world);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(argc != 2)
    {
       if(rank == 0)
            fprintf(stderr,"Usage:Average num_els_per_proc\n");
      exit(1); 
    }

    int num_els_per_proc = atoi(argv[1]);

    double * nums = NULL;
    if(rank == 0)
      nums = createRand(num_els_per_proc * (world));
    
    // make sub arrays
    double * subs = new double[num_els_per_proc];

    // send out all the data to all processes
    MPI_Scatter(nums,num_els_per_proc,MPI_DOUBLE,
            subs,num_els_per_proc,MPI_DOUBLE,
            0,MPI_COMM_WORLD);

    // compute local average
    double avg_sub = computeAvg(subs,num_els_per_proc);

    // gather all partial averages down to root process
    double *sub_avgs = new double[world];

    MPI_Allgather(&avg_sub,1,MPI_DOUBLE,
            sub_avgs,1,MPI_DOUBLE,MPI_COMM_WORLD);

    double avg = computeAvg(sub_avgs,world);
    printf("Process %d,l average = %f\n",rank,avg);

    if(rank == 0)
    {
        delete [] nums;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    delete [] sub_avgs;
    delete [] subs;
    MPI_Finalize();
}


