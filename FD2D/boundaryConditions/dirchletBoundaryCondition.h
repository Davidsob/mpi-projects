//
//  dirchletBoundaryCondition.h
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#ifndef __FD2D__dirchletBoundaryCondition__
#define __FD2D__dirchletBoundaryCondition__

#include <stdio.h>
#include "boundaryCondition.h"

class dirchletBoundaryCondition : public boundaryCondition
{
public:
    
    dirchletBoundaryCondition(double U = 0)
    :boundaryCondition("dirchlet"), U(U){};
    virtual void updateBC(vector<double> &u_ghost);
    virtual void setBC(double &u){ u = this->U;}
private:
    
    double U; // prescribed value
};

#endif /* defined(__FD2D__dirchletBoundaryCondition__) */
