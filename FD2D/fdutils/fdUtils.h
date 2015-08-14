//
//  fdUtils.h
//  FD2D
//
//  Created by Bradley Davidson on 7/31/15.
//
//

#ifndef FD2D_fdUtils_h
#define FD2D_fdUtils_h
#include <math.h>
#include <string>
#include <iostream>

using namespace std;

namespace FDUtils {
    
    typedef enum GRID_DIRECTION{
        
        WEST,
        SOUTH,
        EAST,
        NORTH,
        BOTTOM,
        TOP,
        
    }GRID_DIRECTION;
    
    //  std::string stringForDirection(GRID_DIRECTION direction);
    
    struct Point2d
    {
        int i,j;
    };
    
    struct Stencil
    {
        double O,W,E,S,N,B,T;
        void print(const char * title) const
        {
            printf("%s:{%f, %f, %f, %f, %f, %f, %f}\n",title, O, W, E, S, N, B, T);
        }
    };
    
    int idxFromCoord(int i, int j, int cols);
    
    // means
    double arithemeticMean(double x1, double x2);
    
    double harmonicMean(double x1, double x2);
    
    double geometricMean(double x1, double x2);
}
#endif
