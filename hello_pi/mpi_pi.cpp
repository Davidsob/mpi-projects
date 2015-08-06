#include <mpi.h>
#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

int main(int argc, char ** argv)
{
    if(argc != 2) return 0; 
    long long MAX_PTS = atoll(argv[1]);
   
    // initialize MPI environment
    MPI_Init(NULL,NULL);
    
    double t1 = MPI_Wtime();
    
    // get the number of processes
    int world_size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    // get rank of the process
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    // hot potato
    long long count = 0;
    double x,y;
    
    // perform test
    srand(world_rank+1);
    long long mx_cnt = MAX_PTS/world_size;
    if( world_rank == world_size-1) mx_cnt += (MAX_PTS - world_size*mx_cnt);
    printf("Process %d: testing %lld points...\n",world_rank, mx_cnt);
    for (int i = 0; i < mx_cnt; i++) {
        x = double(rand())/double(RAND_MAX);
        y = double(rand())/double(RAND_MAX);
        if(x*x + y*y <= 1.0) count++;
    }

    if(world_rank != 0)
    {
        MPI_Send(&count,1,MPI_LONG_LONG,0,0,MPI_COMM_WORLD); // send to head node
    }
    
    if(world_rank == 0)
    {
        printf("Process %d has a of %lld. Now summing from other processes...\n",world_rank,count);
        long long count_from;
        for(int i = 1; i < world_size; i ++)
        {
            MPI_Recv(&count_from, 1, MPI_LONG_LONG, i,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            count += count_from;
            printf("Process %d received count %lld from process %d. Count is now %lld\n",world_rank,count_from,i,count);
        }
        
        double pi = 4.0*double(count)/double(MAX_PTS);
        printf("\n\npi calculated to be approx. %f\n",pi);
        double t2 = MPI_Wtime();
        printf("took approximately %f[s] for %d processes and %lld points\n\n\n", t2-t1,world_size, MAX_PTS);
    }

    
    // Finalize the MPI environment
    MPI_Finalize();
}
