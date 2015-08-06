//
//  FDSources.h
//  FD2D
//
//  Created by Bradley Davidson on 7/31/15.
//
//

#ifndef __FD2D__BaseSource__
#define __FD2D__BaseSource__

#include <stdio.h>

namespace FDSource
{
  class PhysicalSource{
  public:

    virtual ~PhysicalSource(){};
    virtual double operator ()(double x, double y, double z , double t) = 0; //functor
  };
}

#endif /* defined(__FD2D__FDSources__) */
