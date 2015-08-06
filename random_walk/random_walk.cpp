#include <mpi.h>
#include <stdlib.h> /* printf, null */
#include <Random.h> /* time */
#include <vector>

using namespace std;

// this is a wrapper around the probe and receive structure
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

void decomposeDomain(int domain_size,
                     int rank,
                     int world_size,
                     int * subdomain_start,
                     int * subdomain_size)
{
  if(world_size > domain_size)
  {
    fprintf(stderr,"Wold size bigger than domain size!\n");
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  *subdomain_size = domain_size/world_size;
  *subdomain_start = *subdomain_size * rank;
  if(rank == world_size-1)
  {
    // the last proccessor gets the remainder
    *subdomain_size += domain_size%world_size;
  }
}

class Walker{
public:
  int location;
  int lifetime;
};


// init the walkers on a process
void initWalkers(int num_walkers_per_proc,
                 int max_walk_size,
                 int subdomain_start,
                 int subdomain_size,
                 vector<Walker> * incoming_walkers)
{
  Walker walker;
  Random rg;
  for(int i = 0; i < num_walkers_per_proc; i++)
  {
    // initialize walker in the middle of the subdomain
    walker.location = subdomain_start;
    walker.lifetime = rg.getRandom(0,max_walk_size);
    incoming_walkers->push_back(walker);
  }
}


void walk(Walker * walker, int subdomain_start, int subdomain_size,
          int domain_size, vector<Walker> * outgoing_walkers){
  while(walker->lifetime > 0)
  {
    if(walker->location == subdomain_start + subdomain_size)
    {
      // take care of case
      // where walker is at end of domain
      // by wrapping it to the beginning
      if(walker->location == domain_size)
      {
        walker->location = 0; // local coord in next domain
      }
      
      outgoing_walkers->push_back(*walker);
      break;
    }else{
      walker->lifetime--;
      walker->location++;
    }
  }
}

void sendOutgoing(vector<Walker> * outgoing_walkers, int rank, int world_size)
{
  int tag =0;
  // clear the outgoing guys
  printf("Process %d sent %lu walkers to process %d\n",
         rank,outgoing_walkers->size(),(rank+1)%world_size);
  
  MPI_Send((void *)outgoing_walkers->data(),
           outgoing_walkers->size()*sizeof(Walker),MPI_BYTE,
           (rank+1)%world_size,tag,MPI_COMM_WORLD);
  

  
  outgoing_walkers->clear();
  
}

void receiveIncoming(vector<Walker> * incoming_walkers, int rank , int world_size)
{
  MPI_Status status;
  int tag = 0;
  int sender_rank = (rank == 0) ? world_size-1 : rank -1;
  
  // use the probe
  MPI_Probe(sender_rank, tag, MPI_COMM_WORLD, &status);
  
  // get the size
  int data_size;
  MPI_Get_count(&status,MPI_BYTE,&data_size);
  
  // set the incoming walkers from raw data
//  vector<Walker> tmp; tmp.reserve(data_size/sizeof(Walker));
  incoming_walkers->resize(data_size/sizeof(Walker));
  // receive the data
  MPI_Recv ((void *)incoming_walkers->data(), data_size, MPI_BYTE,
            sender_rank,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  
//  incoming_walkers->insert(incoming_walkers->begin(), tmp.begin(),tmp.end());
  
  printf("Process %d received %d walkers from process %d\n",
         rank,data_size/sizeof(Walker),sender_rank);
}

int main(int argc, const char *argv[])
{
  // init mpi
  MPI_Init(NULL,NULL);
  
  // domain size
  int domain_size = 20;
  int max_walk_size = 7;
  int walker_per_proc = 1;
  
  if(argc > 1) domain_size = atoi(argv[1]);
  if(argc > 2) max_walk_size = atoi(argv[2]);
  if(argc > 3) walker_per_proc = atoi(argv[3]);
  
  // get processor standing
  int rank = 0, world = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &world);
  
  // get rank of this processor
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // set up the decomposition of my domain
  int subdomain_start,subdomain_size;
  decomposeDomain(domain_size, rank, world, &subdomain_start, &subdomain_size);
  
  // init the walkers in my domain
  vector<Walker> incoming;
  vector<Walker> outgoing;
  initWalkers(walker_per_proc,max_walk_size,
              subdomain_start,subdomain_size, &incoming);
  
  int max_sends_recieves = max_walk_size/(domain_size/world) +1;
  int m = 0;
  while (m < max_sends_recieves) {
    for(int i = 0; i < incoming.size(); i++)
    {
      walk(&incoming[i], subdomain_start,subdomain_size,domain_size,&outgoing);
    }
    
    if(rank %2 == 0)
    {
    // send outgoing
    sendOutgoing(&outgoing,rank,world);
    
    // reveive incoming
    receiveIncoming(&incoming,rank,world);
      
    }else{

      // reveive incoming
      receiveIncoming(&incoming,rank,world);
      // send outgoing
      sendOutgoing(&outgoing,rank,world);
    }
    m++;
  }
  
  // end session
  MPI_Finalize();
}
