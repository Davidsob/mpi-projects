#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

#include <localGrid.h>
#include <fdUtils.h>
#include <boundaryConditions.h>
#include <physicalSources.h>
#include <FDDataManager.h>
#include <mpi.h>

using namespace std;
using namespace FDUtils;
using namespace FDGrid;
using namespace FDSource;

#define MIN(A,B) A < B ? A : B
#define MAX(A,B) A > B ? A : B
#define CLAMP(A,B,C) A > B && A < C

class FDModel{

  public:

    FDModel(int _dim = 0)
    : dim(_dim), CFL(1.0),
    t(0), t_start(0), t_end(0),
    write_every(1){};
    
    virtual ~FDModel(){
    }
    
    virtual void initModel() = 0;
    virtual void applyInitialConditions(MPI_Comm comm) = 0;
    virtual void setBoundaryConditions(MPI_Comm comm){}; // do nothing, but required
    virtual void applyBoundaryConditions(MPI_Comm comm) = 0;
    virtual void updateBoundaryConditions(MPI_Comm comm){};
    virtual void solve(MPI_Comm comm) = 0;
    virtual void advanceSolution(double, MPI_Comm comm) = 0; // advance solution for new time
    virtual void write(MPI_Comm comm){}; // default do nothingi
    virtual void setBoundaryConditions(){} // default do nothing
    virtual double getTimeStep(MPI_Comm comm) = 0;

    // getters

    LocalGrid * getGrid(){return this->grid;}
    double getCFL(){return CFL;}
    double getTime(){return this->t;}
    int getDimension(){return this->dim;}
    size_t static getFileCount(){return FDModel::file_count;}
    
    // setters
    void setDim(int _dim){this->dim = _dim;}
    void setCFL(double _cfl){this->CFL = _cfl;}
    void setTime(double _t){this->t = _t;}
    void setStartTime(double _t0){this->t_start = _t0;}
    void setEndTime(double _tEnd){this->t_end = _tEnd;}
    
    virtual void setGrid(LocalGrid * _grid){
        this->grid = _grid;
    }
    
    virtual void setDataManager(FDDataManager * dmngr){this->data_manager = dmngr;}
    
    void setOuputFilePath(const string &path){this->output_path = path;}
    void setFileName(const string &name){this->base_file_name = name;}
    void setWriteEveryNthStep(size_t n){this->write_every = n;}

    
protected:
    
    int dim;
    double CFL;
    double t;
    double t_start;
    double t_end;
    LocalGrid * grid;
    FDDataManager * data_manager;
    
    // for writing data
    static size_t file_count;
    string output_path;
    string base_file_name;
    size_t write_every;
    
    // utility
    void getDataFromNeighbors(vector<double> &data,
                              vector<vector<double>> &nbrs_data,
                              MPI_Comm comm);
    
    virtual double calculateDivergence(const FDUtils::Stencil &u,
                                        const FDUtils::Stencil &v,
                                        const FDUtils::Stencil &w);
    
    virtual double calculatePartialDerivative(const FDUtils::Stencil &U, int dim);
    
    virtual void updateData(){};
    
    virtual void updateSources();

  private:

};

