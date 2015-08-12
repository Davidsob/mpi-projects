//
//  RadialSource.h
//  FD2D
//
//  Created by Bradley Davidson on 7/31/15.
//
//
#pragma once

#include <stdio.h>
#include "BaseSource.h"
#include <math.h>

namespace FDSource {
    
    // gaussian source functor
    class RadialSource2D : public PhysicalSource{
        
    public:
        RadialSource2D(double xc, double yc, double A, int dim)
        : xC(xc), yC(yc), A(A), dimension(dim) {}
        
        RadialSource2D(){};
        
        virtual double operator ()(double x, double y,double z, double t)
        {
            double theta = atan2(y-yC, x-xC);
            return (dimension == 0) ? -this->A * sin(theta) : this->A * cos(theta);
        }
        
    private:
        double xC,yC;
        double A;
        int dimension;
    };
    
    
}

