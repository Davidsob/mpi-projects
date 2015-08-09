//
//  neumannBoundaryCondition.cpp
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#include "neumannBoundaryCondition.h"

void neumannBoundaryCondition::updateBC(vector<double> &u_ghost)
{
    for(size_t i = 0; i < u_ghost.size(); i++)
    {
        u_ghost[i] = u_ghost[i] - 2.0*this->delX*this->q;
    }
}