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
#include <boundaryCondition.h>

using namespace std;
using namespace FDUtils;

namespace FDGrid
{
    
    class LocalGrid
    {
    public:
        
        /// constructor
        LocalGrid() : ij({0,0}),l_rows(0),l_cols(0), boundary_grid(false){}
        
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
        
        /// return u, uw, us, ue, un // return a stenci
        void stencilPoints(const Point2d &ij,
                           const vector<double> &u,
                           const vector<vector<double>> &nbr_u, FDUtils::Stencil &S, int print = 0) const;
        
        // return boundary points for boundary {
        vector<double> getBoundaryValues(FDUtils::GRID_DIRECTION direction,
                                         const vector<double> &u) const;

        // write grid to file stream
        void write(ofstream &file) const;
        
        // get local 1d idx for grid
        int idxFromCoord(int i, int j) const
        { return i * (l_cols+1) + j;}
        
        // setters
        void setNumberOfRows(int r){this->l_rows = r;};
        void setNumberOfCols(int c){this->l_cols = c;};
        void setColor(int C){ this->color = C;}
        
        // set
        void setNeighbors(const vector<int> &nbrs){
            this->neighbors = nbrs;
            for(int &i : this->neighbors)
            {
                if(i == -1)
                {
                    this->boundary_grid = true;
                    return;
                }
            }
        }
            
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
        
        bool isBoundaryGrid(){return this->boundary_grid;};
        
    private:
        
        struct FDUtils::Point2d ij; // coordinate of local grid in global space
        
        int l_rows,
        l_cols; // cells in  domain
        
        int color; // local process color
        
        bool boundary_grid;
        
        vector<int> neighbors; // neibors list (w,s,e,n)
                
    };
    
}


#endif /* defined(__FD2D__LocalGrid__) */
