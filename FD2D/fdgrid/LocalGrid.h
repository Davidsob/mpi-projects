//
//  LocalGrid.h
//  FD2D
//
//  Created by Bradley Davidson on 8/3/15.
//
//

#ifndef __FD2D__LocalGrid__
#define __FD2D__LocalGrid__

#include <stdio.h>
#include <fstream>
#include <vector>
#include <mpi.h>

#include <fdUtils.h>


using namespace std;
using namespace FDUtils;

namespace FDGrid
{
  
  class LocalGrid
  {
  public:
    
    /// constructor
    LocalGrid() : ij({0,0}),l_rows(0),l_cols(0){}
    
    /// destructor
    ~LocalGrid(){};
    
    /// pass neighbors
    void SendDataToNeighbor(FDUtils::GRID_DIRECTION direction,
                            const vector<double> &u,
                            vector<vector<double>> &nbrs_data,
                            MPI_Comm comm) const;
    /// receive neighbors
    void ReceiveDataFromNeighbor(FDUtils::GRID_DIRECTION direction,
                                 vector<double> &nbr_u,
                                 MPI_Comm comm) const;
    
    /// return u, uw, us, ue, un
    void stencilPoints(const Point2d &ij, double &uW,
                                  double &uS, double &uE, double &uN,
                                  const vector<double> &u,
                                  const vector<vector<double>> &nbr_u, int print = 0) const;
    
    // write grid to file stream
    void write(ofstream &file) const;
    
    // get local 1d idx for grid
    int idxFromCoord(int i, int j) const
    { return i * (l_cols+1) + j;}
    
    // setters
    void setNumberOfRows(int r){this->l_rows = r;};
    void setNumberOfCols(int c){this->l_cols = c;};
    void setColor(int C){ this->color = C;}
    void setNeibors(const vector<int> &nbrs){ this->neighbors = nbrs;}
    void setGlobalCoordinate(int i, int j){ij.i = i, ij.j = j;}
    
    //getters
    int getRows() const
    {return this->l_rows;};
    
    int getCols() const
    {return this->l_cols;};
    
    int getColor() const
    {return this->color;};
    
    Point2d getLocalCoordinate() const
    {
      return this->ij;
    }
    
    const vector<int> & getNeighbors() const
    {return this->neighbors;}
    
  private:
    // public properties    
    struct FDUtils::Point2d ij; // coordinate of local grid in global space
    
    int l_rows,
    l_cols; // cells in  domain
    
    int color; // local process color
    
    vector<int> neighbors; // neibors list (w,s,e,n)
    

    
  };
  
}


#endif /* defined(__FD2D__LocalGrid__) */
