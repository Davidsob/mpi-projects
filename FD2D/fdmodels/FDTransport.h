#pragma once

#include "FDModel.h"

class FDTransport : public FDModel{

  public:
    FDTransport(int dim = 0);  
    virtual ~FDTransport();
    
    virtual void initModel();
    virtual void applyInitialConditions(MPI_Comm comm);
    
    virtual void setBoundaryConditions(vector<boundaryCondition *> _bcs)
    {this->bcs = _bcs;}
    
    virtual void applyBoundaryConditions(MPI_Comm comm);
    virtual void updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm);
    virtual void advanceSolution(double dt,MPI_Comm comm);
    virtual void solve(MPI_Comm comm);
    virtual void write(MPI_Comm comm);
    virtual double getTimeStep();
    
  private:
    
    double calculateAdvection(const FDUtils::Stencil &U,
                              const FDUtils::Stencil &u,
                              const FDUtils::Stencil &v,
                              const FDUtils::Stencil &w);
    
    vector<boundaryCondition *> bcs;
    double calculateMaxU();
    double maxU; // max speed
};
