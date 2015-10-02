#pragma once

#include "FDModel.h"

class FDEuler : public FDModel{
    
public:
    FDEuler(vector<const string> variable_names = {"density", "px", "py","energy"} , int dim = 2);
    virtual ~FDEuler(){};
    
    virtual void initModel();
    virtual void applyInitialConditions(MPI_Comm comm);
    
    /* \brief set a boundary condition for a variable by name.
     * \param[in] name is the name of the variable to apply the boundary condition to
     * \param[in] _bcs is a vector of pointers to boundary condition objects
     *
     * \return \c void
     */
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
    
    virtual void setGamma(double gam){this->gamma = gam;}
    virtual void setGasConstant(double gas){this->R_gas = gas;}
    virtual void setTimeStep(double dt){this->dt = dt;}
    virtual void setArtificialViscosity(vector<double> eps){this->flux_viscocity = eps;}
    void addICName(const string &var, const string &ic_name){
        this->initial_conditions.insert(pair<string, string>(var, ic_name));
    }
    
private:
    
    /// vector of primary variable names
    const vector<const string> primary_variables;
    
    /// vecetor of flux variable names
    vector<string> flux_variables;
    
    /// vector of flux coefficients
    vector<double> flux_viscocity;
    // constants
    double gamma;
    double R_gas;
    double dt;
    
    /// initial condition names
    map<string, string> initial_conditions;
    
    void calculateFlux(const vector<const string> & U);
    void updateData();
    void calculatePressure();
    void calculateTemperature();
    void calculateSpeedOfSound();
    
    //    void calculateGridDivergence();
    //    void calculateMomentum();
    double calculateMaxU();
    
    map<string, vector<boundaryCondition *>> bcs;
};
