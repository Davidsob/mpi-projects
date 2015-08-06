//
//  LocalGrid.cpp
//  FD2D
//
//  Created by Bradley Davidson on 8/3/15.
//
//

#include "LocalGrid.h"
using namespace FDGrid;

#include <string>

#define MIN(A,B) A < B ? A : B
#define MAX(A,B) A > B ? A : B
#define CLAMP(A,B,C) A > B && A < C

std::string stringForDirection(GRID_DIRECTION direction)
{
  std::string dir;
  switch (direction) {
    case FDUtils::WEST: dir = "west"; break;
    case FDUtils::EAST: dir = "east"; break;
    case FDUtils::SOUTH: dir = "south"; break;
    case FDUtils::NORTH: dir = "nort"; break;
      
    default: dir = "none"; break;
  }
  
  return dir;
}

/// pass neighbors
void LocalGrid::SendDataToNeighbor(FDUtils::GRID_DIRECTION direction,
                                   const vector<double> &u,
                                   vector<vector<double>> &nbrs_data,
                                   MPI_Comm comm) const
{
  
  vector<double> tmp;
  int tag = 0;
  
  switch (direction) {
    case FDUtils::WEST:
      for(int i = 0; i <= this->l_rows; i++)
      {
        tmp.push_back(u[this->idxFromCoord(i, 1)]);
      }
      tag = 0;
      break;
      
    case FDUtils::EAST:
      for(int i = 0; i <= this->l_rows; i++)
      {
        tmp.push_back(u[this->idxFromCoord(i, this->l_cols-1)]);
        
      }
      tag = 2;
      break;
      
    case FDUtils::SOUTH:
      for(int i = 0; i <= this->l_cols; i++)
      {
        tmp.push_back(u[this->idxFromCoord(1, i)]);
      }
      tag = 1;
      break;
      
    case FDUtils::NORTH:
      for(int i = 0; i <= this->l_cols; i++)
      {
        tmp.push_back(u[this->idxFromCoord(this->l_rows-1, i)]);
      }
      tag = 3;
      break;
      
    default:
      break;
  }
  
  // which neighbor do we send to
  int send_to = this->neighbors[tag];
  
  // if send to == -1 send to self!
  if(send_to == -1) nbrs_data[tag] = tmp;
  else{
    
    MPI_Send(tmp.data(),tmp.size()*sizeof(double),MPI_BYTE,
             send_to,tag,comm);
  }
  
}
/// receive neighbors
void LocalGrid::ReceiveDataFromNeighbor(FDUtils::GRID_DIRECTION direction,
                                        vector<double> &nbr_u,
                                        MPI_Comm comm) const
{
  int tag = 0;
  int nbr_idx = 0;
  int count = 0;
  switch (direction) {
    case FDUtils::WEST: tag = 2, count = this->l_rows+1, nbr_idx = 0; break;
    case FDUtils::EAST: tag = 0, count = this->l_rows+1, nbr_idx = 2; break;
    case FDUtils::SOUTH: tag = 3, count = this->l_cols+1, nbr_idx = 1; break;
    case FDUtils::NORTH: tag = 1, count = this->l_cols+1, nbr_idx = 3; break;
      
    default:
      break;
  }
  
  int receive_from = this->neighbors[nbr_idx];
  if (receive_from == -1) return;
  if(receive_from == -1) MPI_Comm_rank(comm, &receive_from);
  
  // resize and receive
  //  int rank;
  //  MPI_Comm_rank(comm, &rank);
  //  printf("Proc %d: receiving data %s, from proc %d with tag %d\n",rank,
  //         stringForDirection(direction).c_str(),
  //         receive_from,tag);
  //  printf("my nbrs are: {%d, %d, %d, %d}\n",
  //         neighbors[0],neighbors[1],neighbors[2],neighbors[3]);
  nbr_u.clear();
  nbr_u.resize(count);
  MPI_Recv(nbr_u.data(),count*sizeof(double),MPI_BYTE,
           receive_from,tag,
           comm,MPI_STATUS_IGNORE);
  
}

void LocalGrid::stencilPoints(const Point2d &ij, double &uW,
                   double &uS, double &uE, double &uN,
                              const vector<double> &u,
                   const vector<vector<double>> &nbr_u, int print) const
{
  if(CLAMP(ij.i, 0, this->l_rows) && CLAMP(ij.j, 0, this->l_cols))
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
  }else if (CLAMP(ij.i, 0, this->l_rows) && ij.j == 0) // west boundary
  {
    uW = nbr_u[0][ij.i];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
    if(print) printf("{%d, %d} = west boundary\n",ij.i,ij.j);
  }else if (ij.i == 0 && CLAMP(ij.j, 0, this->l_cols)) // south boundary
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = nbr_u[1][ij.j];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
    if(print) printf("{%d, %d} = south boundary\n",ij.i,ij.j);
  }else if (CLAMP(ij.i, 0, this->l_rows) && ij.j == this->l_cols) // east boundary
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = nbr_u[2][ij.i];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
    if(print) printf("{%d, %d} = east boundary\n",ij.i,ij.j);
  }else if (ij.i == this->l_rows && CLAMP(ij.j, 0, this->l_cols)) // north boundary
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = nbr_u[3][ij.j];
    if(print) printf("{%d, %d} = north boundary\n",ij.i,ij.j);
  }else if (ij.i == 0 && ij.j == 0) // southwest corner
  {
    uW = nbr_u[0][ij.i];
    uS = nbr_u[1][ij.j];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
    if(print) printf("{%d, %d} = sw corner\n",ij.i,ij.j);
  }else if (ij.i == 0 && ij.j == this->l_cols) // southeast corner
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = nbr_u[1][ij.j];
    uE = nbr_u[2][ij.i];
    uN = u[idxFromCoord(ij.i+1, ij.j)];
    if(print) printf("{%d, %d} = se corner\n",ij.i,ij.j);
  }else if (ij.i == this->l_rows && ij.j == this->l_cols) // northeast corner
  {
    uW = u[idxFromCoord(ij.i, ij.j-1)];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = nbr_u[2][ij.i];
    uN = nbr_u[3][ij.j];
    if(print) printf("{%d, %d} = ne corner\n",ij.i,ij.j);
  }else if (ij.i == this->l_rows && ij.j == 0) // northwest corner
  {
    uW = nbr_u[0][ij.i];
    uS = u[idxFromCoord(ij.i-1, ij.j)];
    uE = u[idxFromCoord(ij.i, ij.j+1)];
    uN = nbr_u[3][ij.j];
    if(print) printf("{%d, %d} = nw corner\n",ij.i,ij.j);
  }else{
    printf("ERROR!!!!!! point {%d, %d} not catergorized!!!!!\n"
           "{rows,cols} = {%d, %d}\n", ij.i, ij.j, l_rows, l_cols);
  }
  
  
}

// write grid to file stream
void LocalGrid::write(ofstream &file) const
{
  file << this->ij.i << "\t";
  file << this->ij.j << "\t";
  file << this->l_rows << "\t";
  file << this->l_cols << "\t";
  file << this->color << "\n";
}