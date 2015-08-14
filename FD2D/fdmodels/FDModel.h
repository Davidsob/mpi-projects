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
    hx(1), hy(1), hz(1), grid (new LocalGrid), write_every(1){};
    
    virtual ~FDModel(){
        delete this->grid;
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
    virtual double getTimeStep() = 0;
    
    // getters
    bool hasData(const string &name)
    {
        return this->data.find(name) != this->data.end();
    }
    
    bool hasSource(const string &name)
    {
        return this->sources.find(name) != this->sources.end();
    }
    
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
    
    void setGridSpacing(double _h, int dim = 0){
        switch (dim) {
            case 1: this->hy = _h; break;
            case 2: this->hz = _h; break;
            default: this->hx = _h; break;
        }
    }
    
    void setData(const string &name, vector<double> &x)
    {
        this->data.insert(std::pair<string, vector<double>>(name,x));
        
        // automatically add shared data
        vector<vector<double>> shared_zeros(4,{0});
        this->shared_data.insert(std::pair<string, vector<vector<double>> >(name,shared_zeros));
        
        // set stencil automatically
        Stencil s;
        this->stencils.insert(std::pair<string, FDUtils::Stencil >(name,s));
    }
    
    void setSource(const string &name, PhysicalSource * S, bool initData = false)
    {
        this->sources.insert(std::pair<string, PhysicalSource *>(name,S));
        if(initData)
        {
            int local_els = (grid->getRows()+1)*(grid->getCols()+1);
            vector<double> tmp(local_els,0);
            this->setData(name, tmp);
        }
    }
    
    void setGrid(LocalGrid * _grid){this->grid = _grid;}
    void setOuputFilePath(const string &path){this->output_path = path;}
    void setFileName(const string &name){this->base_file_name = name;}
    void setWriteEveryNthStep(size_t n){this->write_every = n;}

    
  protected:
    int dim;
    double CFL;
    double t;
    double t_start;
    double t_end;
    double hx,hy,hz; // uniform grid spacing
    
    LocalGrid * grid;
    map<string, vector<double>> data;
    map<string, PhysicalSource *> sources;
    map<string, vector<vector<double>> > shared_data; //neighboring data
    map<string, FDUtils::Stencil> stencils; //stencil structs
    
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
    
    vector<double> &getData(const string &name)
    {
        return this->data[name];
    }
    
    vector<vector<double>> &getSharedData(const string &name)
    {
         return this->shared_data[name];
    }
    
    FDUtils::Stencil &getStencil(const string &name)
    {
         return this->stencils[name];
    }
    
    PhysicalSource * getSource(const string &name)
    {
         return this->sources[name];
    }
    
    virtual void updateData(){};
    
  private:

};

