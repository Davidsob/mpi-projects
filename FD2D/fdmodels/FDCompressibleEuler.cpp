#include "FDCompressibleEuler.h"
#include "FDHeatTransfer.h"
#include "FDTransport.h"
#include <math.h>

static  vector< string> flux0{"E0", "E1", "E2", "E3"};
static  vector< string> flux1{"F0", "F1", "F2", "F3"};
static  vector< string> flux2{"G0", "G1", "G2", "G3"};

FDEuler::FDEuler(vector<const string> variable_names, int dim)
: FDModel(dim), primary_variables(variable_names), gamma(7.0/5.0), R_gas(1.0), dt(1.0e-3)
{
    this->flux_variables.insert(this->flux_variables.end(),flux0.begin(),flux0.end());
    
    if (dim > 1)
        this->flux_variables.insert(this->flux_variables.end(),flux1.begin(),flux1.end());
    
    if(dim > 2)
        this->flux_variables.insert(this->flux_variables.end(),flux2.begin(),flux2.end());
    
}


void FDEuler::initModel()
{
    // set local storage
    int local_els = this->grid->getNumberOfGridPoints();
    vector<double> zeros(local_els,0);
    for(auto var : primary_variables)
    {
        this->data_manager->setData(var, zeros);
        this->data_manager->setData(var + "_*", zeros);
    }
    
    // init Fluxes
    for(auto var : flux_variables)
    {
        this->data_manager->setData(var, zeros);
    }
    
    this->data_manager->setData("pressure", zeros);
    this->data_manager->setData("temperature", zeros);
    this->data_manager->setData("c", zeros);
}


void FDEuler::applyBoundaryConditions(MPI_Comm comm)
{

    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    // get vector of neighboring nodes
    vector<int> nbr = this->grid->getNeighbors();
    
    // set dirchlet
    for(auto bc : this->bcs)
    {
        vector<double> &Ui = this->data_manager->getData(bc.first);
        if(nbr[0] == -1)
            for(i = 0; i <= local_rows; i++) bc.second[0]->setBC(Ui[this->grid->idxFromCoord(i, 0)]);
        if(nbr[1] == -1)
            for(i = 0; i <= local_cols; i++) bc.second[1]->setBC(Ui[this->grid->idxFromCoord(0, i)]);
        if(nbr[2] == -1)
            for(i = 0; i <= local_rows; i++) bc.second[2]->setBC(Ui[this->grid->idxFromCoord(i, local_cols)]);
        if(nbr[3] == -1)
            for(i = 0; i <= local_cols; i++) bc.second[3]->setBC(Ui[this->grid->idxFromCoord(local_rows, i)]);
        
    }
    MPI_Barrier(comm);
}

void FDEuler::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
    const vector<int> &nbrs = this->grid->getNeighbors();
    
    for(auto bc : this->bcs)
    {
        vector<vector<double>> & nbr_data = this->data_manager->getSharedData(bc.first);
        vector<double> &Ui = this->data_manager->getData(bc.first);
        for(size_t i = 0; i < nbrs.size(); i++)
        {
            if(nbrs[i] == -1)
            {
                if(bc.second[i]->getName() == "convectiveCooling")
                {
                    vector<double> boundary_val = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i), Ui);
                    dynamic_cast<convectiveCooling *>(bc.second[i])->updateBC(nbr_data[i], boundary_val);
                }else{
                    bc.second[i]->updateBC(nbr_data[i]);
                }
            }
        }
    }
}

void FDEuler::applyInitialConditions(MPI_Comm comm)
{
    // want to get left and right values from neihbors
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    for(auto icpair : initial_conditions)
    {
        // Get data reference from manager
        vector<double> &U = this->data_manager->getData(icpair.first);
        
        // get ic for this data
        PhysicalSource * ic = this->data_manager->getSource(icpair.second);
        
        // apply IC function for x and y (currently only support 2d)
        double x, y;
        for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++) {
            x = this->grid->getX(idx);
            y = this->grid->getY(idx);
            U[idx] = ic->operator()(x, y, 0, this->t_start);
        }
    }
    
    
    MPI_Barrier(comm);
}

void FDEuler::solve(MPI_Comm comm)
{
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    // 0. update initial state of sources
//    this->updateSources();
    size_t step_no = 0;
    // write initial solution
    this->updateData(); // calculate secondary variables
    this->write(comm);
    
    MPI_Barrier(comm);
    
    // begin solve
    double dt = 0;
    while (this->t <= this->t_end) {
        
        //0. set time
        dt = this->getTimeStep(comm);
        this->setTime(this->t + dt);
        
        //1. update all sources (some may be time dependent)
        this->updateSources();
        
        step_no++;
//        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        // 2. Advnace solution one time step using MacCormack Scheme
        this->advanceSolution(dt,comm);

        // 3. write data to file
        MPI_Barrier(comm);
        if((step_no % write_every) == 0 || this->t >= this->t_end)
        {
            if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
            this->updateData(); // calculate secondary variables
            this->write(comm);
        }
        
    }
    
}

void FDEuler::advanceSolution(double dt,MPI_Comm comm)
{
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    // Advance the time solution Via MacCormack's method
    // share data
    for (auto p : this->data_manager->availableData())
        FDModel::getDataFromNeighbors(this->data_manager->getData(p),
                                      this->data_manager->getSharedData(p), comm);
    // hold up
    MPI_Barrier(comm);
    
    // update boundary conditions
    for(auto var : primary_variables)
        this->updateBoundaryConditions(this->data_manager->getSharedData(var),comm);
    
    size_t idx;
    Point2d ij{0,0};
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    // go nuts
    vector<double> &rho_s = this->data_manager->getData(primary_variables[0]+ "_*");
    vector<double> &px_s = this->data_manager->getData(primary_variables[1]+ "_*");
    vector<double> &py_s = this->data_manager->getData(primary_variables[2]+ "_*");
    vector<double> &e_s = this->data_manager->getData(primary_variables[3]+ "_*");
    
    //  vars
    vector<double> &rho = this->data_manager->getData(primary_variables[0]);
    vector<double> &px = this->data_manager->getData(primary_variables[1]);
    vector<double> &py = this->data_manager->getData(primary_variables[2]);
    vector<double> &e = this->data_manager->getData(primary_variables[3]);
    
    double Fx = 0, Fy = 0;
    double dtdx = dt / this->grid->get_hx();
    double dtdy = dt / this->grid->get_hy();

    // MacCormack Step 1
    // calculate flux
    this->calculateFlux({"density", "px", "py", "energy"});
    
    for(int i = 0; i <= local_rows; i++)
    {
        for(int j = 0; j <= local_cols; j++)
        {
            ij.i = i; ij.j = j;
            idx = this->grid->idxFromCoord(i , j);
            
            // set all stencils
            for (auto p : this->data_manager->availableData())
            {
                this->grid->stencilPoints(ij, this->data_manager->getData(p),
                                          this->data_manager->getSharedData(p),
                                          this->data_manager->getStencil(p));
            }
            
            // update u
            /// diffusion term
            Stencil E0 = this->data_manager->getStencil("E0");
            Stencil E1 = this->data_manager->getStencil("E1");
            Stencil E2 = this->data_manager->getStencil("E2");
            Stencil E3 = this->data_manager->getStencil("E3");
            Stencil F0 = this->data_manager->getStencil("F0");
            Stencil F1 = this->data_manager->getStencil("F1");
            Stencil F2 = this->data_manager->getStencil("F2");
            Stencil F3 = this->data_manager->getStencil("F3");
            
            // gravitational body force
            if(this->data_manager->hasSource("gx"))
                Fx = this->data_manager->getStencil("gx").O * rho[idx] ;
            
            // gravitational body force
            if(this->data_manager->hasSource("gy"))
                Fy = this->data_manager->getStencil("gy").O * rho[idx] ;
            
            
            // MacCormack Step 1
            rho_s[idx] = rho[idx] - dtdx * (E0.E - E0.O) - dtdy * (F0.N - F0.O);
            px_s[idx] = px[idx] - dtdx * (E1.E - E1.O) - dtdy * (F1.N - F1.O);
            py_s[idx] = py[idx] - dtdx * (E2.E - E2.O) - dtdy * (F2.N - F2.O);
            e_s[idx] = e[idx] - dtdx * (E3.E - E3.O) - dtdy * (F3.N - F3.O);
            
        }
    }
    
    // MacCormack Step 2
    // calculate flux
    this->calculateFlux({"density_*", "px_*", "py_*", "energy_*"});
    
    for(int i = 0; i <= local_rows; i++)
    {
        for(int j = 0; j <= local_cols; j++)
        {
            ij.i = i; ij.j = j;
            idx = this->grid->idxFromCoord(i , j);
            
            // set all stencils
            for (auto p : this->data_manager->availableData())
            {
                this->grid->stencilPoints(ij, this->data_manager->getData(p),
                                          this->data_manager->getSharedData(p),
                                          this->data_manager->getStencil(p));
            }
            
            // update u
            /// diffusion term
            Stencil E0 = this->data_manager->getStencil("E0");
            Stencil E1 = this->data_manager->getStencil("E1");
            Stencil E2 = this->data_manager->getStencil("E2");
            Stencil E3 = this->data_manager->getStencil("E3");
            Stencil F0 = this->data_manager->getStencil("F0");
            Stencil F1 = this->data_manager->getStencil("F1");
            Stencil F2 = this->data_manager->getStencil("F2");
            Stencil F3 = this->data_manager->getStencil("F3");
            
            // gravitational body force
            if(this->data_manager->hasSource("gx"))
                Fx = this->data_manager->getStencil("gx").O * rho[idx] ;
            
            // gravitational body force
            if(this->data_manager->hasSource("gy"))
                Fy = this->data_manager->getStencil("gy").O * rho[idx] ;
            
            
            // MacCormack Step 1
            rho[idx] = 0.5*((rho[idx] + rho_s[idx]) - dtdx * (E0.O - E0.W) - dtdy * (F0.O - F0.S));
            px[idx] = 0.5*((px[idx] + px[idx]) - dtdx * (E1.O - E1.W) - dtdy * (F1.O - F1.S));
            py[idx] = 0.5*((py[idx] + py[idx]) - dtdx * (E2.O - E2.W) - dtdy * (F2.O - F2.S));
            e[idx] = 0.5*((e[idx] + e[idx]) - dtdx * (E3.O - E3.W) - dtdy * (F3.O - F3.S));
            
        }
    }
}

void FDEuler::write(MPI_Comm comm)
{
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    // write to file
    char title[50];
    static int file_no = -1;
    
    file_no++;
    sprintf(title, "/%s.%03i.proc_%d",this->base_file_name.c_str(),file_no,rank);
    ofstream file;
    file.open(this->output_path + string(title));
    
    //    printf("Data to write:\n");
    //    for(auto s : this->data_manager->availableData())
    //        printf("data soure: %s\n", s.c_str());
    vector<double> &rho = this->data_manager->getData(primary_variables[0]);
    vector<double> &px = this->data_manager->getData(primary_variables[1]);
    vector<double> &py = this->data_manager->getData(primary_variables[2]);
    vector<double> &e = this->data_manager->getData(primary_variables[3]);
    vector<double> &T = this->data_manager->getData("temperature");
    vector<double> &p = this->data_manager->getData("pressure");
    vector<double> &c = this->data_manager->getData("c");
    
    
    double x,y;
    for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++)
    {
        x = this->grid->getX(idx);
        y = this->grid->getY(idx);
        file << x << "\t" << y;
        file << "\t" << rho[idx] <<"\t" << px[idx] << "\t" << py[idx] << "\t" << p[idx];
        file << "\t" << e[idx] << "\t" << T[idx] <<"\t" << c[idx] << "\n";
    }
    
    file.close();
    if(rank == 0)
    {
        FDModel::file_count++;
        sprintf(title, "/time.dat");
        if(this->t == this->t_start)
            file.open(this->output_path + string(title));
        else
            file.open(this->output_path + string(title),ios::app);
        
        file << this->t <<"\n";
        file.close();
    }
}


// A method for calculating secondary variables (fluxes, pressure, temperature, etc)
void FDEuler::calculateFlux(const vector<const string> &U)
{
    
    vector<double> &rho = this->data_manager->getData(U[0]);
    vector<double> &px = this->data_manager->getData(U[1]);
    vector<double> &py = this->data_manager->getData(U[2]);
    vector<double> &e = this->data_manager->getData(U[3]);
    
    vector<double> &E0 = this->data_manager->getData("E0");
    vector<double> &E1 = this->data_manager->getData("E1");
    vector<double> &E2 = this->data_manager->getData("E2");
    vector<double> &E3 = this->data_manager->getData("E3");
    
    vector<double> &F0 = this->data_manager->getData("F0");
    vector<double> &F1 = this->data_manager->getData("F1");
    vector<double> &F2 = this->data_manager->getData("F2");
    vector<double> &F3 = this->data_manager->getData("F3");
    
    double fact, p;
    for(size_t idx = 0; idx < rho.size(); idx++){
        p= (this->gamma - 1.0) * (e[idx] - 0.5 * (pow(px[idx], 2) + pow(py[idx], 2))/rho[idx]);
        
        // E
        fact = px[idx]/rho[idx];
        E0[idx] = fact * rho[idx];
        E1[idx] = fact * px[idx] + p;
        E2[idx] = fact * py[idx];
        E3[idx] = fact * (e[idx] + p);
        
        // F
        fact = py[idx]/rho[idx];
        F0[idx] = fact * rho[idx];
        F1[idx] = fact * px[idx];
        F2[idx] = fact * py[idx] + p;
        F3[idx] = fact * (e[idx] + p);
    }
}

void FDEuler::updateData()
{
    // update data for export
    this->calculatePressure();
    this->calculateTemperature();
    this->calculateSpeedOfSound();
}


void FDEuler::calculatePressure()
{
    vector<double> &rho = this->data_manager->getData("density");
    vector<double> &px = this->data_manager->getData("px");
    vector<double> &py = this->data_manager->getData("py");
    vector<double> &e = this->data_manager->getData("energy");
    vector<double> &p = this->data_manager->getData("pressure");
    
    for(size_t idx = 0; idx < rho.size(); idx++){
        p[idx] = (this->gamma - 1.0) * ( e[idx] - 0.5 * (pow(px[idx], 2) + pow(py[idx], 2))/rho[idx]);
    }
}

void FDEuler::calculateTemperature()
{
    vector<double> &rho = this->data_manager->getData("density");
    vector<double> &p = this->data_manager->getData("pressure");
    vector<double> &T = this->data_manager->getData("temperature");

    for(size_t idx = 0; idx < rho.size(); idx++){
        T[idx] = this->R_gas * this->gamma * p[idx]/(rho[idx] * (this->gamma - 1.0));
    }
}

void FDEuler::calculateSpeedOfSound()
{
    vector<double> &c = this->data_manager->getData("c");
    vector<double> &T = this->data_manager->getData("temperature");
    
    for(size_t idx = 0; idx < T.size(); idx++){
        c[idx] = sqrt(this->R_gas * this->gamma * T[idx]);
    }
}

double FDEuler::getTimeStep(MPI_Comm comm)
{
//    FDModel * m = getModel("mass transport");
//    return dynamic_cast<FDTransport *>(m)->getTimeStep(comm);
    return this->dt;
}

