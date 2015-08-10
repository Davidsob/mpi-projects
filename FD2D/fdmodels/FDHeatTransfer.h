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
    virtual void solve(MPI_Comm comm);
    virtual void write(MPI_Comm comm);
    virtual double getTimeStep();
    
    void setDensity(PhysicalSource * _rho){this->rho = _rho;}
    void setHeatCapacity(PhysicalSource * _Cp){this->Cp = _Cp;}
    void setThermalConductivity(PhysicalSource * _K){this->K = _K;}
    void set_u(PhysicalSource * _u){this->u = _u;}
    void set_v(PhysicalSource * _v){this->v = _v;}
    void set_w(PhysicalSource * _w){this->w = _w;}
    void set_p(PhysicalSource * _p){this->p = _p;}
    
  private:
    
    double calculateDiffusion(const FDUtils::Stencil &T, const FDUtils:: Stencil &K);
    double calculateAdvection(const FDUtils::Stencil &T,
                              const FDUtils::Stencil &u,
                              const FDUtils::Stencil &v,
                              const FDUtils::Stencil &w);
    
    double calculateDeformationEnergy(const FDUtils::Stencil &T,const FDUtils::Stencil &p);
    


    PhysicalSource * rho; // density
    PhysicalSource * Cp; // constant pressure heat capacity
    PhysicalSource * K; // thermal conductivity
    PhysicalSource * u; // u velocity component
    PhysicalSource * v; // v velocity component
    PhysicalSource * w; // w velocity component
    PhysicalSource *p; // pressure
    
    vector<boundaryCondition *> bcs;

};
