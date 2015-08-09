//
//  dirchletBoundaryCondition.cpp
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#include "dirchletBoundaryCondition.h"

void dirchletBoundaryCondition::updateBC(vector<double> &u_ghost)
{
    for(size_t i = 0; i < u_ghost.size(); i++)
    {
        u_ghost[i] = 2.0*this->U - u_ghost[i];
    }
}
