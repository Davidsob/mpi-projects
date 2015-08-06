//
//  GaussianSource.h
//  FD2D
//
//  Created by Bradley Davidson on 7/31/15.
//
//

#ifndef __FD2D__GaussianSource__
#define __FD2D__GaussianSource__

#include <stdio.h>
#include "BaseSource.h"
#include <math.h>

namespace FDSource {
  
  // gaussian source functor
  class GaussianSource : public PhysicalSource{
    
  public:
    GaussianSource(double xc, double yc, double sigx, double sigy, double A)
    : xc(xc), yc(yc), sigx(sigx), sigy(sigy), A(A){}
    
    ~GaussianSource(){};
    
    virtual double operator ()(double x, double y,double z, double t)
    {
      double s = A*exp(-0.5*(pow((x-xc)/sigx, 2.0) + pow((y-yc)/sigy, 2.0)));
      return s;
    }
    
  private:
    double xc,yc,sigx,sigy,A;
  };
  
}

#endif /* defined(__FD2D__GaussianSource__) */
