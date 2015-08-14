#pragma once

#include "FDModel.h"


class FDCompressibleFlow : public FDModel{
    
public:
    FDCompressibleFlow(int dim = 0);
    virtual ~FDCompressibleFlow();
    
    virtual void initModel();
    virtual void applyInitialConditions(MPI_Comm comm);
    
    virtual void setBoundaryConditions(const string &name, vector<boundaryCondition *> _bcs)
    {
        this->bcs.insert(pair<string, vector<boundaryCondition *>>(name,_bcs));
    }
    
    virtual void applyBoundaryConditions(MPI_Comm comm);
    virtual void updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm);
    virtual void advanceSolution(double dt,MPI_Comm comm);
    virtual void solve(MPI_Comm comm);
    virtual void write(MPI_Comm comm);
    virtual double getTimeStep();
    
    FDModel * getModel(const string & name)
    {
        return this->couplings[name];
    }
    
private:
    
    map<string, FDModel *> couplings;
    
    void addModel(const string &name, FDModel * m)
    {
        this->couplings.insert(pair<string, FDModel *>(name,m));
    }
    double calculateDiffusion(const FDUtils::Stencil &T, const FDUtils:: Stencil &K);
    double calculateAdvection(const FDUtils::Stencil &T,
                              const FDUtils::Stencil &u,
                              const FDUtils::Stencil &v,
                              const FDUtils::Stencil &w,
                              const FDUtils::Stencil &div);
    
    double calculateDeformationEnergy(const FDUtils::Stencil &T,const FDUtils::Stencil &p);
    
    void updateData();
    void calculatePressure();
    void calculateGridDivergence();
    void calculateMomentum();
    
    map<string, vector<boundaryCondition *>> bcs;
};
