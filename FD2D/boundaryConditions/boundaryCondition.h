//
//  boundaryCondition.h
//  FD2D
//
//  Created by Bradley Davidson on 8/7/15.
//
//

#ifndef FD2D_boundaryCondition_h
#define FD2D_boundaryCondition_h

#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
using namespace std;

class boundaryCondition
{
public:
    
    boundaryCondition(string kind = "base") : kind(kind){};
    virtual ~boundaryCondition(){};
    virtual void updateBC(vector<double> &u_ghost)
    {
        std::runtime_error("base boundary condition updateBC not implemented");
    } // do nothing
    
    virtual void setBC(double &u){}; // by default do nothing
    virtual const string &getName(){return this->kind;}
    
protected:
    string kind;
};
#endif
