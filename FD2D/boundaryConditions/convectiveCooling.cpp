//
//  convectiveCooling.cpp
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#include "convectiveCooling.h"

void convectiveCooling::updateBC(vector<double> &u_ghost,const vector<double> &u_boundary)
{
    double q = 0;
    for(size_t i = 0; i < u_ghost.size(); i++)
    {
        q = this->h * (this->u_infinity - u_boundary[i]);
        u_ghost[i] = u_ghost[i] - 2.0*this->delX*q;
    }
}