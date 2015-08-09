#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <mpi.h>

#include <communicatorMap.h>
#include <LocalGrid.h>
#include <GaussianSource.h>
#include <fdUtils.h>

using namespace std;
using namespace FDSource;
using namespace FDGrid;
using namespace FDUtils;

#define MIN(A,B) A < B ? A : B
#define MAX(A,B) A > B ? A : B
#define CLAMP(A,B,C) A > B && A < C

static int filecount = 0;
static const string saveto = "/Users/LNLB/mpi-projects/FD2D/";

// set initial conditions for problem
void setIC(vector<double> & um, vector<double> &u0,
           double hx, double hy,double C,
           PhysicalSource * H,
           LocalGrid &grid, MPI_Comm comm);

// file write time loop
void timeLoop(vector<double> &um, vector<double> &u, vector<double> &up,
              double hx, double hy, double t_end, double C,
              PhysicalSource * H,
              LocalGrid &grid,MPI_Comm comm);

// write grid partition to file
void writeGridPartition(const vector<LocalGrid> & parts)
{
  ofstream file;
  file.open(saveto + string("grids/graph.dat"));
  for(const LocalGrid & p : parts) p.write(file);
  file.close();
}

// write data to file
void writeToFile2D(const vector<double> & u,
                   double hx, double hy, double t,
                   const LocalGrid & grid,
                   MPI_Comm comm);

void getWorldInfo(MPI_Comm comm, int &r, int &ws)
{
  MPI_Comm_rank(comm,&r);
  MPI_Comm_size(comm,&ws);
}


int main(int argc, char ** argv)
{
  MPI_Init(NULL,NULL);
  
  double timer = 0.0;
  timer -= MPI_Wtime();
  
  int rank, nproc;
  getWorldInfo(MPI_COMM_WORLD,rank,nproc);
  
  int divx = 15;
  int divy = 15;
  int Nx = divx+1;
  int Ny = divy+1;
  
  
  // Partition Global Grid
  LocalGrid my_grid;
  vector<LocalGrid> mesh;
  vector<int> all_neighbors;
  vector<LocalGrid>grids;
  if(rank == 0)
  {
    CommunicatorMap map(nproc);
    
    for(int i = 0; i < nproc; i++)
    {
      LocalGrid lg;
      map.setLocalGrid(Nx, Ny, i, &lg);
      grids.push_back(lg);
    }
    
    // scatter nbrs
    all_neighbors = map.stackNbrs();
    
    // print the grid partition
    writeGridPartition(grids);
  }
  
  //  // scatter local grids
  MPI_Scatter(grids.data(),sizeof(LocalGrid),MPI_BYTE,
              &my_grid, sizeof(LocalGrid),MPI_BYTE,
              0,MPI_COMM_WORLD);
  
  //  // scatter neighbor info to local grids
  vector<int> tmp(4,0);
  MPI_Scatter(all_neighbors.data(),4*sizeof(int),MPI_BYTE,
              tmp.data(), 4*sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);
  
  my_grid.setNeighbors(tmp);
  MPI_Barrier(MPI_COMM_WORLD);
  
  // print results on local processors
  //  printf("proc %d:\n [%d, %d]x[%d, %d]\ncolor: %i\n",rank,my_grid.idx0,my_grid.idx0 +my_grid.ino-1,
  //         my_grid.jdx0, my_grid.jdx0+my_grid.jno-1,my_grid.color);
  //
  //  if(!my_grid.neighbors.empty())
  //    printf("neighbors: %d, %d %d,%d\n\n",
  //           my_grid.neighbors[0], my_grid.neighbors[1],
  //           my_grid.neighbors[2],my_grid.neighbors[3]);
  
  // set local storage
  int local_els = (my_grid.getRows()+1)*(my_grid.getCols()+1);
  vector<double> um(local_els,0);
  vector<double> u(local_els,0);
  vector<double> up(local_els,0);

  // get source
  GaussianSource * H = new GaussianSource(0.5,0.5,0.1,0.1,0.1);
  
  // Courant number
  double C;
  C = 0.05;
  
  double t_end = 10.0;
  double hx = 1.0/(Nx-1.0) , hy = 1.0/(Ny - 1.0);
  
  // set IC
  setIC(um,u,hx,hy,C,dynamic_cast<PhysicalSource *>(H), my_grid, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // write to file
  timeLoop(um,u,up,hx, hy, t_end, C,
           dynamic_cast<PhysicalSource *>(H),
           my_grid, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // write file count
  if(rank == 0)
  {
    double dt =  MIN(hx, hy)*C;
    ofstream file;
    file.open(saveto + "solutions/model_info.dat");
    file << filecount << "\n";
    file << nproc << "\n";
    file << C << "\n";
    file << dt <<"\n";
    file.close();
    
    timer += MPI_Wtime();
    printf("Simulation complete: time/proc = %f[s]\n",timer/nproc);
    printf("Simulation complete....\n");
  }
  
  MPI_Finalize();
}

void getDataFromNeighbors(vector<vector<double>> &nbrs_data,
                          const vector<double> &u,
                          const LocalGrid &grid, MPI_Comm comm)
{
  // pass receive depending on color
  
  if(grid.getColor() == 0)
  {
    // send
    grid.SendDataToNeighbor(WEST,u,nbrs_data,comm);
    grid.SendDataToNeighbor(SOUTH,u,nbrs_data,comm);
    grid.SendDataToNeighbor(EAST,u,nbrs_data,comm);
    grid.SendDataToNeighbor(NORTH,u,nbrs_data,comm);
    // receive
    grid.ReceiveDataFromNeighbor(WEST,nbrs_data[0],comm);
    grid.ReceiveDataFromNeighbor(SOUTH,nbrs_data[1],comm);
    grid.ReceiveDataFromNeighbor(EAST,nbrs_data[2],comm);
    grid.ReceiveDataFromNeighbor(NORTH,nbrs_data[3],comm);
  }else if(grid.getColor() == 1){
    
    // receive
    grid.ReceiveDataFromNeighbor(WEST,nbrs_data[0],comm);
    grid.ReceiveDataFromNeighbor(SOUTH,nbrs_data[1],comm);
    grid.ReceiveDataFromNeighbor(EAST,nbrs_data[2],comm);
    grid.ReceiveDataFromNeighbor(NORTH,nbrs_data[3],comm);
    
    // send
    grid.SendDataToNeighbor(WEST,u,nbrs_data,comm);
    grid.SendDataToNeighbor(SOUTH,u,nbrs_data,comm);
    grid.SendDataToNeighbor(EAST,u,nbrs_data,comm);
    grid.SendDataToNeighbor(NORTH,u,nbrs_data,comm);
  }
}


void setIC(vector<double> & um, vector<double> &u0,
           double hx, double hy,double C, PhysicalSource * H,
           LocalGrid &grid, MPI_Comm comm)
{
  
  // want to get left and right values from neihbors
  int rank , nproc;
  getWorldInfo(MPI_COMM_WORLD,rank,nproc);
  
  // get local informaion
  int local_rows = grid.getRows();
  int local_cols = grid.getCols();
  Point2d local_ij = grid.getLocalCoordinate();
  
  size_t i, j, idx;
  double dt = MIN(hx, hy)*C, rx = dt/hx, ry = dt/hy;
  double x,y;
  
  // 1 set initial displacement
  vector<double> Hij(u0.size(),0);
  for (i = 0; i <= local_rows; i++) {
    for(j = 0; j <= local_cols; j++)
    {
      idx = grid.idxFromCoord(i,j);
      x = (local_ij.j + j)*hx;
      y = (local_ij.i + i)*hy;
      Hij[idx] = H->operator()(x,y, 0.0,0.0);
      u0[idx] = Hij[idx];
    }
  }
  
  MPI_Barrier(comm);
  // 2.get neighboring data points
  vector<vector<double>> neighbors_data(4,{0});
  getDataFromNeighbors(neighbors_data,u0,grid,comm);
  
  //  // set um
  double rxsq = rx*rx;
  double rysq = ry*ry;
  double d2x, d2y;
  
  Point2d ij{0,0};
  // note hij = u0ij
  double uij = 0, uW = 0, uE = 0, uS = 0, uN = 0;
  for(int i = 0; i <= local_rows; i++)
  {
    for(j = 0; j <= local_cols; j++)
    {
      ij.i = i;
      ij.j = j;
      idx = grid.idxFromCoord(i, j);
      uij = u0[idx];
      
      grid.stencilPoints(ij, uW, uS, uE, uN, u0, neighbors_data);
      
      d2x =  FDUtils::arithemeticMean(uE,uij)*(uE - uij) -
      FDUtils::arithemeticMean(uW, uij)*(uij - uW);
      
      d2y = FDUtils::arithemeticMean(uN,uij)*(uN - uij) -
      FDUtils::arithemeticMean(uS, uij)*(uij - uS);
      
      um[idx] = uij + 0.5 * (rxsq * d2x + rysq * d2y);
    }
  }
  
  //    um[0] = 0; um[n-1] = 0;
}

void timeLoop(vector<double> & um, vector<double> &u, vector<double> &up,
              double hx, double hy, double t_end, double C,
              PhysicalSource * H,
              LocalGrid &grid,MPI_Comm comm)
{
  
  int rank , nproc;
  getWorldInfo(MPI_COMM_WORLD,rank,nproc);
  
  size_t i, j, idx, step_no = 0;
  double dt = MIN(hx, hy)*C, rx = dt/hx, ry = dt/hy, t = 0;
  
  writeToFile2D(u,hx,hy,0,grid,MPI_COMM_WORLD);
  
  // get local informaion
  int local_rows = grid.getRows();
  int local_cols = grid.getCols();
  Point2d local_ij = grid.getLocalCoordinate();
  
  // get
  double x,y;
  vector<double> Hij(u.size(),0);
  for (i = 0; i <= local_rows; i++) {
    for(j = 0; j <= local_cols; j++)
    {
      idx = grid.idxFromCoord(i,j);
      x = (local_ij.j + j)*hx;
      y = (local_ij.i + i)*hy;
//      Hij[idx] = 1.0 -  H->operator()(x,y, 0.0,t);
      Hij[idx] = 1.0;
    }
  }
  
  MPI_Barrier(comm);

  double uij = 0, uW = 0, uE = 0, uS = 0, uN = 0;
  double hij = 0, hW = 0, hE = 0, hS = 0, hN = 0;
  double d2x, d2y;

  double rxsq = rx*rx;
  double rysq = ry*ry;
  vector<vector<double>> neighbors_data(4,{0});
  vector<vector<double>> neighbors_source(4,{0});
  Point2d ij{0,0};
  while (t <= t_end) {
    t += dt;
    step_no++;
    // share data
    // 2.get neighboring data points
    getDataFromNeighbors(neighbors_data,u,grid,comm);
    MPI_Barrier(comm);
    
    getDataFromNeighbors(neighbors_source,Hij,grid,comm);
    MPI_Barrier(comm);
    
    for(int i = 0; i <= local_rows; i++)
    {
      for(j = 0; j <= local_cols; j++)
      {
        ij.i = i;
        ij.j = j;
        idx = grid.idxFromCoord(i , j);
        uij = u[idx];
        
        grid.stencilPoints(ij, uW, uS, uE, uN, u, neighbors_data);
        grid.stencilPoints(ij, hW, hS, hE, hN, Hij, neighbors_source);
        
        d2x =  FDUtils::arithemeticMean(hE,hij)*(uE - uij) -
        FDUtils::arithemeticMean(hW, hij)*(uij - uW);
        
        d2y = FDUtils::arithemeticMean(hN,hij)*(uN - uij) -
        FDUtils::arithemeticMean(hS, hij)*(uij - uS);
        
        up[idx] = 2.0 * uij - um[idx] + (rxsq * d2x + rysq * d2y);
      
      }
    }

    um = u;
    u = up;
    MPI_Barrier(comm);
    writeToFile2D(u,hx,hy,t,grid,MPI_COMM_WORLD);
  }
  
}

void writeToFile2D(const vector<double> & u,
                   double hx, double hy, double t,
                   const  LocalGrid & grid,
                   MPI_Comm comm)
{
  // want to get left and right values from neihbors
  int rank , nproc;
  getWorldInfo(MPI_COMM_WORLD,rank,nproc);
  
  // write to file
  char title[50];
  static int f = -1;
  
  f++;
  sprintf(title, "solutions/u.dat.%03i.proc_%d",f,rank);
  ofstream file;
  file.open(saveto + string(title));
  
  
  double x,y;
  int idx;
  int i,j;
  Point2d local_ij = grid.getLocalCoordinate();
  int rows = grid.getRows();
  int cols = grid.getCols();
  
  for (i = 0; i <= rows; i++) {
    for(j = 0; j <= cols; j++)
    {
      idx = grid.idxFromCoord(i, j);
      x = (local_ij.j + j)*hx;
      y = (local_ij.i + i)*hy;
      file << x << "\t" << y << "\t" << u[idx] <<"\n";
    }
  }
  
  file.close();
  
  if(rank == 0) filecount++;
}
