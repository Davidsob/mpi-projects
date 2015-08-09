//
//  convectiveCooling.h
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#ifndef __FD2D__convectiveCooling__
#define __FD2D__convectiveCooling__

#include <stdio.h>
#include "boundaryCondition.h"

class convectiveCooling : public boundaryCondition
{
public:
    convectiveCooling(double h = 0, double u_infinity = 0, double delX = 0)
    :boundaryCondition("convectiveCooling"), h(h),u_infinity(u_infinity), delX(delX){};
    
    virtual void updateBC(vector<double> &u_ghost, const vector<double> &u_boundary);
private:
    double h; // convective coefficient
    double u_infinity; // prescribed flux
    double delX;
};


#endif /* defined(__FD2D__convectiveCooling__) */
