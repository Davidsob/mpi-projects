#pragma once

#include "FDModel.h"

class FDHeatTransfer : public FDModel{

  public:
    FDHeatTransfer(int dim = 0);  
    virtual ~FDHeatTransfer();
    
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
    
    double calculateDiffusion(const FDUtils::Stencil &T, const FDUtils:: Stencil &K);
    double calculateAdvection(const FDUtils::Stencil &T,
                              const FDUtils::Stencil &u,
                              const FDUtils::Stencil &v,
                              const FDUtils::Stencil &w);
    
    double calculateDeformationEnergy(const FDUtils::Stencil &T,const FDUtils::Stencil &p);
    
    vector<boundaryCondition *> bcs;

};
