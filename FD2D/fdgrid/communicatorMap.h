#ifndef COMMUNICATORMAP
#define COMMUNICATORMAP

#include <stdio.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "LocalGrid.h"
#include <fdUtils.h>

using namespace std;
using namespace FDUtils;
// grid stuff

namespace FDGrid
{
  class CommunicatorMap
  {
  public:
    // default constructor
    CommunicatorMap(int proc = 0);
    
    // destructor
    ~CommunicatorMap();
    
    // partition the grid
    void partitionGrid(int N);
    
    // get idx from coords
    int idxFromCoord(int i, int j){ return i * this->cols + j;}
    
    vector<int> stackNbrs(); // return a stacked array of all neighbors

    // getters
    /// get ij coordinate of a process in the process grid space
    Point2d getCoordinateOfProcess(int p){
      return this->proc_to_grid[p];
    }
    
    /// get color of process p
    int getColorOfProcess(int p){return colors[p];}
    
    /// get array of neghbors for process p
    vector<int> getNbr(int p){return nbrs[p];}
    
    /// set up local grid for some x-y space
    // completes decomposition of domain
    void setLocalGrid(int Nx, int Ny, int rank,LocalGrid * localGrid);
    
  private:
    
    int rows,
        cols; // grows and cols of communicator grid
    
    Point2d * proc_to_grid; // map of proc to grid
    vector<vector<int>> nbrs; // nbrs for all processors
    vector<int> colors; // colors for all processors
    
    void setNbr(int i, int j, vector<int> &nbr); // set a processors neighbor
  };
  
}


#endif
