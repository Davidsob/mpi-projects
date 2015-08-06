#include "communicatorMap.h"

using namespace FDGrid;

CommunicatorMap::CommunicatorMap(int proc)
{
  
  partitionGrid(proc);
  
  // hash addresses/ set nbrs / set colors
  size_t idx;
  vector<int> tmp(4,-1);
  this->nbrs = vector<vector<int>>(proc,tmp);
  this->colors = vector<int>(proc,-1);
  proc_to_grid = new struct Point2d[proc];
  
  for(int i = 0; i < this->rows; i ++)
    for(int j = 0; j < this->cols; j ++){
      idx = this->idxFromCoord(i,j);
      proc_to_grid[idx] = Point2d{i,j};
      setNbr(i, j, nbrs[idx]);
    }
  
}

CommunicatorMap::~CommunicatorMap(){
  if(proc_to_grid)
    delete [] proc_to_grid;
}

void CommunicatorMap::partitionGrid(int N)
{
  this->rows = sqrt(N);
  this->cols = N/this->rows;
  
  while (N%this->cols != 0 && this->rows > 0) {
    this->rows--;
    this->cols = N/this->rows;
  }
  
  if(this->rows == 0)
  {
    this->rows = 1;
    this->cols = N;
  }
}

void CommunicatorMap::setNbr(int i, int j, vector<int> &nbr)
{
  int idx;
  int my_proc = idxFromCoord(i,j);
  if(my_proc == 0) colors[0] = 0;
  
  if(j > 0){
    idx = idxFromCoord(i, j-1);
    nbr[0] = idx;
    if(colors[idx] == -1) colors[idx] = (colors[(my_proc)]+1)%2;
  }
  
  if(i > 0)
  {
    idx = idxFromCoord(i-1, j);
    nbr[1] = idx;
    if(colors[idx] == -1) colors[idx] = (colors[(my_proc)]+1)%2;
  }
  
  if(j < this->cols-1){
    idx = idxFromCoord(i, j+1);
    nbr[2] = idx;
    if(colors[idx] == -1) colors[idx] = (colors[(my_proc)]+1)%2;
  }
  
  if(i < this->rows-1){
    idx = idxFromCoord(i+1, j);
    nbr[3] = idx;
    if(colors[idx] == -1) colors[idx] =(colors[(my_proc)]+1)%2;
  }
}

vector<int> CommunicatorMap::stackNbrs(){
  vector<int> tmp;
  for (vector<int> &v : this->nbrs) {
    tmp.insert(tmp.end(), v.begin(), v.end());
  }
  return tmp;
}

void CommunicatorMap::setLocalGrid(int Nx, int Ny,int rank, LocalGrid * localGrid)
{
  //  if(rows > Nx || cols > Ny)
  //  {
  //    fprintf(stderr,"Wold size bigger than domain size!\n");
  ////    MPI_Abort(MPI_COMM_WORLD,1);
  //  }

  // get coordinate in proc grid
  int local_rows, local_cols, local_i, local_j;
  Point2d ij = this->getCoordinateOfProcess(rank);
    
  local_rows = (Ny-1)/this->rows;
  local_i = local_rows * ij.i;
  
  if(ij.i == this->rows-1)
    local_rows += (Ny-1)%this->rows;
  
  local_cols = (Nx-1)/this->cols;
  local_j = local_cols * ij.j;
  
  if(ij.j == this->cols-1)
    local_cols += (Nx-1)%this->cols;
  
  // set color
  localGrid->setGlobalCoordinate(local_i, local_j);
  localGrid->setNumberOfRows(local_rows);
  localGrid->setNumberOfCols(local_cols);
  localGrid->setColor(this->getColorOfProcess(rank));
}

