 //
//  main.cpp
//  FD1D
//
//  Created by Bradley Davidson on 7/26/15.
//  Copyright (c) 2015 PhunPhactory. All rights reserved.
//

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <mpi.h>

using namespace std;
static int filecount = 0;

static const string saveto = "/Users/davidson/MPI-Tutorials/FD1D/solutions/";

void setIC(vector<double> &u0, vector<double> & um,
           int idxa, double h, double C, MPI_Comm comm);

void timeLoop(vector<double> &up,vector<double> &u, vector<double> & um,
              double h, double t_end, double C,MPI_Comm comm);

void writeToFile(const vector<double> & u, double t, MPI_Comm comm);


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

int main(int argc, const char * argv[]) {
  
  // init the mpi world
  MPI_Init(NULL,NULL);
  
  double timer = 0.0;
  timer -= MPI_Wtime();
  // find standing on processor
  int rank , world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  // Model params
  int div;
  double C, t_end;
  div =argc > 1 ? atoi(argv[1]) : 20;
  C = argc > 2 ? atof(argv[2]) : 0.3;
  t_end = argc > 3 ? atof(argv[3]) : 1.0;
  int nx = div+1;
  double h = 1.0/div;
  if(rank ==0) // print model
    printf("\n\nMODEL INFO:\ncells: %i\nnx = %d\nhx = %f\nCFL: %f\nt_end: %f[s]\n\n",
           div,nx,h,C, t_end);

  // get partition of mesh for proc
  int idxa, localLen;
  decomposeDomain(nx, rank, world_size, &idxa, &localLen);
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // set asside some memory
  vector<double> up(localLen,0);
  vector<double> u(localLen,0);
  vector<double> um(localLen,0);
  
  //  // set IC
  setIC(u,um,idxa,h,C,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // sync before starting time loop
  
  //
  //  // solve
  timeLoop(up,u,um,h,t_end,C,MPI_COMM_WORLD);
  
  // write file count and model info
  if(rank == 0)
  {
    ofstream file;
    file.open(saveto + "file_count.txt");
    file << filecount;
    file <<"\n" << C;
    file << "\n" << nx-1 <<"\n";
    file.close();
  }
  // cout << "...Simulation complete." << endl;
  // printf("Process %d, idx = %d, idxb = %d\n",rank,idxa,idxb);
  
  MPI_Barrier(MPI_COMM_WORLD);
  timer += MPI_Wtime();
  if(rank == 0) printf("Simulation complete: time/proc = %f[s]\n",timer/world_size);
  MPI_Finalize();
}

void getNbrNode(const vector<double> &u,int rank,int world_size, MPI_Comm comm,
                double &left, double &right)
{
  // get neihbors info
  if(rank % 2 == 0)
  {
    if(rank != world_size-1)
      MPI_Send(const_cast<double *>(&u[u.size()-1]),1,MPI_DOUBLE,rank+1,0,comm);
    if(rank != 0)
      MPI_Send(const_cast<double *>(&u[0]),1,MPI_DOUBLE,rank-1,0,comm);
    
    if(rank != world_size-1)
      MPI_Recv(&right, 1, MPI_DOUBLE,rank+1,0,comm,MPI_STATUS_IGNORE);
    if(rank != 0)
      MPI_Recv(&left, 1, MPI_DOUBLE,rank-1,0,comm,MPI_STATUS_IGNORE);
    
  }else{
    
    if(rank != world_size-1)
      MPI_Recv(&right, 1, MPI_DOUBLE,rank+1,0,comm,MPI_STATUS_IGNORE);
    if(rank != 0)
      MPI_Recv(&left, 1, MPI_DOUBLE,rank-1,0,comm,MPI_STATUS_IGNORE);
    
    if(rank != world_size-1)
      MPI_Send(const_cast<double *>(&u[u.size()-1]),1,MPI_DOUBLE,rank+1,0,comm);
    if(rank != 0)
      MPI_Send(const_cast<double *>(&u[0]),1,MPI_DOUBLE,rank-1,0,comm);
  }
  
}

void setIC(vector<double> &u0, vector<double> & um,
           int idxa, double h, double C, MPI_Comm comm)
{
  
  // want to get left and right values from neihbors
  int rank , world_size;
  MPI_Comm_size(comm, &world_size);
  MPI_Comm_rank(comm, &rank);
  
  // set initial displacement
  double x;
  double umx = 0.05, csq = C*C;
  // set initial displacement
  x = idxa*h;
  
  for(double &ui : u0)
  {
    if(x < 0.7) ui = umx/0.7 * x;
    else ui = umx/0.3*(1.0-x);
    x += h;
  }
  
  MPI_Barrier(comm); // sync before passing neighbor info
  
  // get neihbors info
  double left = 0, right = 0;
  getNbrNode(u0,rank,world_size,comm,left,right);
  
  // set um
  size_t local_len = u0.size();
  for(size_t i = 0; i < local_len ; i++)
  {
    if(i == 0)
      um[i] = u0[i] + 0.5*csq*(u0[i+1]-2*u0[i]+left);
    else if(i == local_len-1)
      um[i] = u0[i] + 0.5*csq*(right-2*u0[i]+u0[i-1]);
    else // interior nodes
      um[i] = u0[i] + 0.5*csq*(u0[i+1]-2*u0[i]+left);
  }
  
  MPI_Barrier(comm); // sync before passing neighbor info
  // apply boundary conditions at the ends
  if(rank == 0) um[0] = left;
  if(rank == world_size-1) um[local_len-1] = 0;
  
}

void timeLoop(vector<double> &up,vector<double> &u, vector<double> & um,
              double h,double t_end, double C, MPI_Comm comm)
{
  
  // want to get left and right values from neihbors
  int rank , world_size;
  MPI_Comm_size(comm, &world_size);
  MPI_Comm_rank(comm, &rank);
  
  size_t i, step_no = 0, loc_nx = u.size();
  double dt = C*h, t = 0, Csq = C*C;
  
  // practice gather and write!
  writeToFile(u,0,MPI_COMM_WORLD);
  
  double left = 0, right = 0;
  while (t <= t_end) {
    t += dt;
    step_no++;
    
    // get neihbors info
    getNbrNode(u,rank,world_size,comm,left,right);
    // update proc cells
    
    for(i = 0; i < loc_nx ; i++)
    {
      
      if(i == 0)
        up[i] = 2*u[i] - um[i] + Csq*(u[i+1]-2*u[i]+left);
      else if(i == loc_nx-1)
        up[i] = 2*u[i] - um[i] + Csq*(right-2*u[i]+u[i-1]);
      else // interior nodes
        up[i] = 2*u[i] - um[i] + Csq*(u[i+1]-2*u[i]+u[i-1]);
      
    }
    
    // apply boundary conditions at the ends
    if(rank == 0) up[0] = left;
    if(rank == world_size-1) up[loc_nx-1] = 0;
    
    um = u;
    u = up;
    
    MPI_Barrier(comm); // sync before write
    writeToFile(u, t, comm);
  }
}

void writeToFile(const vector<double> & u, double t, MPI_Comm comm)
{
  // want to get left and right values from neihbors
  int rank , world_size;
  MPI_Comm_size(comm, &world_size);
  MPI_Comm_rank(comm, &rank);
  
  // get total size
  size_t * local_nx = NULL;
  if(rank == 0)
  {
    local_nx = new size_t[world_size];
  }
  
  size_t nxl = u.size();
  MPI_Gather(&nxl,1,MPI_UNSIGNED_LONG,
             local_nx,1, MPI_UNSIGNED_LONG,0, MPI_COMM_WORLD);
  MPI_Barrier(comm); // sync before write
  
  double * uall = NULL;
  size_t nx = 0;
  if(rank == 0)
  {
    for(int k = 0; k < world_size; k++) nx += local_nx[k];
    uall = new double[nx+1];
  }
  
  MPI_Gather(const_cast<double *>(&u[0]),nxl,MPI_DOUBLE,
             uall,nxl, MPI_DOUBLE,0, MPI_COMM_WORLD);
  
  if (rank == 0) {
    double h = 1.0/(nx-1);
    char title[50];
    static int i = -1;
    
    i++;
    sprintf(title, "u.dat.%03i",i);
    ofstream file;
    file.open(saveto + string(title));
    double x = 0;
    for(size_t j = 0; j < nx; j++)
    {
      file << x << "\t" << uall[j] << "\n";
      x += h;
    }
    
    file.close();
    filecount++;
    delete [] uall;
    delete [] local_nx;
  }
  
}
