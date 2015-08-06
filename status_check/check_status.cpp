#include <mpi.h>
#include <stdlib.h> /* printf, null */
#include <Random.h> /* time */

void MPI_ProbeRecv(void ** data,
                    MPI_Datatype datatype,
                    int source,
                    int tag,
                    MPI_Comm comm,
                    MPI_Status * status)
{
    MPI_Probe(source,tag,comm,status);
    int data_size = 0;
    MPI_Get_count(status, datatype, &data_size);
    void * tmp = (void *)malloc(data_size*sizeof(datatype));
    MPI_Recv(tmp,data_size,datatype,source,tag,comm,status);
    *data = tmp;
}


int main(int argc, const char *argv[])
{
    // init mpi
    MPI_Init(NULL,NULL);
    
    // get processor standing
    int rank = 0, world = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &world);
    
    // get rank of this processor
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // send some random numbers
    int number_amount;
    if(rank == 0)
    {
        const int MAX_NUMBERS = 10;
        int numbers[MAX_NUMBERS];

        Random rg;
        number_amount = rg.getRandom(0,MAX_NUMBERS-1);
        
        // make numbers
        for(int i = 0; i < number_amount; i++)
          numbers[i] = rg.getRandom(-10,10);
        // send
        MPI_Send(numbers,number_amount,MPI_INT,1,0,MPI_COMM_WORLD);
        printf("0 sends %d numbers to 1\n",number_amount);
        
        for(int i = 0; i < number_amount; i++)
          printf("%d, ", numbers[i]);
        printf("\n");      

    }else if(rank == 1)
    {
        MPI_Status status;
        
        // probe for message and size
       // MPI_Probe(0,0,MPI_COMM_WORLD,&status);
        
        // check status out to find out how many numbers were actually sent
        //MPI_Get_count(&status, MPI_INT, &number_amount);
        
        //int * numbers = new int[number_amount];
        // recieve from 0
        //MPI_Recv(numbers,number_amount,MPI_INT,0,0,MPI_COMM_WORLD, &status);
        void * numbers = NULL;
        MPI_ProbeRecv(&numbers,MPI_INT,0,0,MPI_COMM_WORLD,&status);

        int number_amount = 0;
        MPI_Get_count(&status,MPI_INT,&number_amount);
        printf("1 received %d numbers from 0. Message source = %d, "
                "tag = %d\n",number_amount, status.MPI_SOURCE, status.MPI_TAG);
        
        // print numbers
        printf("Try to print numbers...\n");
        int * tmp = (int *)numbers;
        for(int i = 0; i < number_amount; i++)
          printf("%d, ", tmp[i]);
        printf("\n"); 
    
        free(numbers);     
    }
    
    // end session
    MPI_Finalize();
}
