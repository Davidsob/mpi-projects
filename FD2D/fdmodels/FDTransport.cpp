#include "FDTransport.h"

FDTransport::FDTransport(int dim)
: FDModel(dim)
{
    bcs = {nullptr,nullptr,nullptr,nullptr};
    maxU = 1;
}

FDTransport::~FDTransport()
{
    
}

void FDTransport::initModel()
{
    // set local storage
    int local_els = (grid->getRows()+1)*(grid->getCols()+1);
    vector<double> zeros(local_els,0);
    
    this->setData("U", zeros);
    this->setData("Up", zeros);
}

void FDTransport::applyBoundaryConditions(MPI_Comm comm)
{
    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    vector<int> nbr = this->grid->getNeighbors();
    vector<double> &U = this->getData("U");
    // set dirchlet
    if(nbr[0] == -1)
        for(i = 0; i <= local_rows; i++) this->bcs[0]->setBC(U[this->grid->idxFromCoord(i, 0)]);
    if(nbr[1] == -1)
        for(i = 0; i <= local_cols; i++) this->bcs[1]->setBC(U[this->grid->idxFromCoord(0, i)]);
    if(nbr[2] == -1)
        for(i = 0; i <= local_rows; i++) this->bcs[2]->setBC(U[this->grid->idxFromCoord(i, local_cols)]);
    if(nbr[3] == -1)
        for(i = 0; i <= local_cols; i++) this->bcs[3]->setBC(U[this->grid->idxFromCoord(local_rows, i)]);
    
    MPI_Barrier(comm);
}

void FDTransport::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
    if(!this->grid->isBoundaryGrid()) return;
    
    const vector<int> &nbrs = this->grid->getNeighbors();
    vector<double> &U = this->getData("U");
    
    for(size_t i = 0; i < nbrs.size(); i++)
    {
        if(nbrs[i] == -1)
        {
            if(this->bcs[i]->getName() == "convectiveCooling")
            {
                vector<double> boundary_u = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i),U);
                dynamic_cast<convectiveCooling *>(bcs[i])->updateBC(nbrs_data[i], boundary_u);
            }else{
                this->bcs[i]->updateBC(nbrs_data[i]);
            }
        }
    }
}

void FDTransport::applyInitialConditions(MPI_Comm comm)
{
    if (!hasSource("initial condition")) return;
    // want to get left and right values from neihbors
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    // get local informaion
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    Point2d local_ij = this->grid->getLocalCoordinate();
    
    size_t i, j, idx;
    double x,y;
    vector<double> &U = this->getData("U");
    PhysicalSource *ic = getSource("initial condition");
    // 1 set initial Temperature at 0
    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            U[idx] = ic->operator()(x, y, 0, this->t_start);;
        }
    }
    
    MPI_Barrier(comm);
}

double FDTransport::calculateAdvection(const FDUtils::Stencil &U,
                                       const FDUtils::Stencil &u,
                                       const FDUtils::Stencil &v,
                                       const FDUtils::Stencil &w)
{
    double advection = 0;
    if(this->dim >= 0)
    {
        if(u.O >= 0)
            advection += u.O * (U.O - U.W)/this->hx;
        else
            advection += u.O * (U.E - U.O)/this->hx;
    }
    
    if(this->dim >= 1)
    {
        if(v.O >= 0)
            advection += v.O * (U.O - U.S)/this->hy;
        else
            advection += v.O * (U.N - U.O)/this->hy;
    }
    
    if(this->dim >= 2)
    {
        if(w.O >= 0)
            advection += w.O * (U.O - U.B)/this->hz;
        else
            advection += w.O * (U.T - U.O)/this->hz;
    }
    
    advection += U.O * this->calculateDivergence(u,v,w);
    return -1.0 * advection;
}

void FDTransport::solve(MPI_Comm comm)
{
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    size_t i, j, idx, step_no = 0;
    
    this->write(comm);
    
    // get local informaion
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    Point2d local_ij = this->grid->getLocalCoordinate();
    
    // get
    double x,y;
    
    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            
            for(auto p : this->sources)
            {
                if(hasData(p.first))
                    getData(p.first)[idx] = p.second->operator()(x, y, 0, this->t);
            }
            
        }
    }
    
    MPI_Barrier(comm);
    
    double dt = this->getTimeStep();
    while (this->t <= this->t_end) {
        
        // advance time step
        this->t += dt;
        step_no++;
        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        this->advanceSolution(dt, comm);
        
        MPI_Barrier(comm);
        if((step_no % write_every) == 0 || this->t >= this->t_end)
            this->write(comm);
    }
    
}

void FDTransport::advanceSolution(double dt, MPI_Comm comm)
{
    // share data
    for (auto p : this->data)
        FDModel::getDataFromNeighbors(getData(p.first), getSharedData(p.first), comm);
    
    // update boundary conditions
    this->updateBoundaryConditions(getSharedData("U"),comm);
    MPI_Barrier(comm);
    
    vector<double> &Up = getData("Up");
    vector<double> &U = getData("U");
    double source = 0;
    size_t idx;
    Point2d ij{0,0};
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    for(int i = 0; i <= local_rows; i++)
    {
        for(int j = 0; j <= local_cols; j++)
        {
            ij.i = i; ij.j = j;
            
            for (auto p : this->stencils)
                this->grid->stencilPoints(ij, getData(p.first), getSharedData(p.first), getStencil(p.first));
            
            /// advection term
            if(hasSource("u") || hasSource("v")|| hasSource("w"))
                source = this->calculateAdvection(getStencil("U"),
                                                  getStencil("u"),
                                                  getStencil("v"),
                                                  getStencil("w"));
            
            idx = this->grid->idxFromCoord(i , j);
            Up[idx] =  U[idx] + dt * source;
            
        }
    }
    
    U = Up;
}

void FDTransport::write(MPI_Comm comm)
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
    
    
    double x,y;
    int idx;
    int i,j;
    Point2d local_ij = this->grid->getLocalCoordinate();
    int rows = this->grid->getRows();
    int cols = this->grid->getCols();
    
    vector<double> &U = this->getData("U");
    for (i = 0; i <= rows; i++) {
        for(j = 0; j <= cols; j++)
        {
            idx = this->grid->idxFromCoord(i, j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            file << x << "\t" << y << "\t" << U[idx] <<"\n";
        }
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

double FDTransport::getTimeStep()
{
    return MIN(this->hx/maxU, this->hy/maxU)*this->CFL;
}

