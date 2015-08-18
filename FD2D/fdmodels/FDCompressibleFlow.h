#pragma once

#include "FDModel.h"

class FDCompressibleFlow : public FDModel{
    
public:
    FDCompressibleFlow(vector<const string> variable_names = {"u", "v", "w"} , int dim = 0);
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
    virtual double getTimeStep(MPI_Comm comm);
    
    FDModel * getModel(const string & name)
    {
        return this->couplings[name];
    }
        
    void setCFL(double _cfl){
        this->CFL = _cfl;
        for(auto m : couplings) m.second->setCFL(_cfl);
    }
    
    void setTime(double _t){
        this->t = _t;
        for(auto m : couplings) m.second->setTime(_t);
    }
    
    virtual void setGrid(LocalGrid * _grid);
    
    virtual void setDataManager(FDDataManager * dmngr);
    
    void addICName(const string &var, const string &ic_name){
        this->initial_conditions.insert(pair<string, string>(var, ic_name));
    }

private:
    
    const vector<const string> primary_variables;
    map<string, string> initial_conditions;
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
