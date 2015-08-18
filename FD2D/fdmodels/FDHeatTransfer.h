#pragma once

#include "FDModel.h"

class FDHeatTransfer : public FDModel{
  
public:
  FDHeatTransfer(string variable_name = "temperature", int dim = 0);
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
  virtual double getTimeStep(MPI_Comm comm);
  void setICName(const string &ic_name){this->initial_condition = ic_name;}
  void setCoupled(bool coupled){ this->coupled = coupled;}
  
private:
  
  string primary_variable;
  string initial_condition;
  bool coupled;
  double calculateDiffusion(const FDUtils::Stencil &T, const FDUtils:: Stencil &K);
  double calculateAdvection(const FDUtils::Stencil &T,
                            const FDUtils::Stencil &u,
                            const FDUtils::Stencil &v,
                            const FDUtils::Stencil &w);
  
  double calculateViscousHeating();
  
  vector<boundaryCondition *> bcs;
};
