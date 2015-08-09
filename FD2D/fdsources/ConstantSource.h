//
//  ConstantSource.h
//  FD2D
//
//  Created by Bradley Davidson on 7/31/15.
//
//

#ifndef __FD2D__ConstantSource__
#define __FD2D__ConstantSource__

#include <stdio.h>
#include "BaseSource.h"
#include <math.h>

namespace FDSource {
  
  // gaussian source functor
  class ConstantSource : public PhysicalSource{
    
  public:
    ConstantSource(double A)
    : A(A){}
    
    ~ConstantSource(){};
    
    virtual double operator ()(double x, double y,double z, double t)
    {
      return this->A;
    }
    
  private:
    double A;
  };
  
}

#endif /* defined(__FD2D__ConstantSource__) */
