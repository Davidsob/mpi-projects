//
//  neumannBoundaryCondition.h
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#ifndef __FD2D__neumannBoundaryCondition__
#define __FD2D__neumannBoundaryCondition__

#include <stdio.h>
#include "boundaryCondition.h"

class neumannBoundaryCondition : public boundaryCondition
{
public:
    neumannBoundaryCondition(double q = 0, double delX = 0)
    : boundaryCondition("neumann"), q(q), delX(delX){};
    virtual void updateBC(vector<double> &u_ghost);
private:
    
    double q; // prescribed flux
    double delX;
};


#endif /* defined(__FD2D__neumannBoundaryCondition__) */
