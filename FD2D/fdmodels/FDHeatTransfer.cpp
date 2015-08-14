#include "FDHeatTransfer.h"

FDHeatTransfer::FDHeatTransfer(int dim)
: FDModel(dim)
{
    bcs = {nullptr,nullptr,nullptr,nullptr};
}

FDHeatTransfer::~FDHeatTransfer()
{

}

void FDHeatTransfer::initModel()
{
    // set local storage
    int local_els = (grid->getRows()+1)*(grid->getCols()+1);
    vector<double> zeros(local_els,0);
    this->setData("T", zeros);
    this->setData("Tp", zeros);
}

void FDHeatTransfer::applyBoundaryConditions(MPI_Comm comm)
{
    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    vector<int> nbr = this->grid->getNeighbors();
    vector<double> &U = this->getData("T");
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

void FDHeatTransfer::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
    if(!this->grid->isBoundaryGrid()) return;
    
    const vector<int> &nbrs = this->grid->getNeighbors();
    
    for(size_t i = 0; i < nbrs.size(); i++)
    {
        if(nbrs[i] == -1)
        {
            if(this->bcs[i]->getName() == "convectiveCooling")
            {
                vector<double> boundary_u = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i), getData("T"));
                dynamic_cast<convectiveCooling *>(bcs[i])->updateBC(nbrs_data[i], boundary_u);
            }else{
                this->bcs[i]->updateBC(nbrs_data[i]);
            }
        }
    }
}

void FDHeatTransfer::applyInitialConditions(MPI_Comm comm)
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
    
    // 1 set initial Temperature at 0
    vector<double> &U = getData("T");
    PhysicalSource * ic = getSource("initial condition");
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

double FDHeatTransfer::calculateDiffusion(const FDUtils::Stencil &T, const FDUtils::Stencil &K)
{
    double diff = 0.0;
    if(this->dim >= 1)
    {
        double rxsq = 1.0/this->hx/this->hx;
        double d2x =  FDUtils::arithemeticMean(K.E,K.O)*(T.E - T.O) -
        FDUtils::arithemeticMean(K.W, K.O)*(T.O - T.W);
        diff += rxsq * d2x;
    }
    
    if(this->dim >= 2)
    {
        double rysq = 1.0/this->hy/this->hy;
        double d2y = FDUtils::arithemeticMean(K.N,K.O)*(T.N - T.O) -
        FDUtils::arithemeticMean(K.S, K.O)*(T.O - T.S);
        diff += rysq * d2y;
    }
    
    if(this->dim == 3)
    {
        double rzsq = 1.0/this->hz/this->hz;
        double d2z = FDUtils::arithemeticMean(K.T,K.O)*(T.T - T.O) -
        FDUtils::arithemeticMean(K.B, K.O)*(T.O - T.B);
        diff += rzsq * d2z;
    }
    
    return diff;
}


double FDHeatTransfer::calculateAdvection(const FDUtils::Stencil &T,
                                          const FDUtils::Stencil &u,
                                          const FDUtils::Stencil &v,
                                          const FDUtils::Stencil &w)
{
    double advection = 0;
    if(this->dim >= 0)
    {
        if(u.O >= 0)
            advection += u.O * (T.O - T.W)/this->hx;
        else
            advection += u.O * (T.E - T.O)/this->hx;
    }
    
    if(this->dim >= 1)
    {
        if(v.O >= 0)
            advection += v.O * (T.O - T.S)/this->hy;
        else
            advection += v.O * (T.N - T.O)/this->hy;
    }
    
    if(this->dim >= 2)
    {
        if(w.O >= 0)
            advection += w.O * (T.O - T.B)/this->hz;
        else
            advection += w.O * (T.T - T.O)/this->hz;
    }
    
    advection += T.O * this->calculateDivergence(u,v,w);
    return -1.0 * advection;
}


double FDHeatTransfer::calculateDeformationEnergy(const FDUtils::Stencil &p,const FDUtils::Stencil &u)
{
    double deformation = 0;
    return deformation;
}

void FDHeatTransfer::solve(MPI_Comm comm)
{
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    size_t i, j, idx, step_no = 0;
    // write initial solution
    this->write(comm);
    
    // get local informaion
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    Point2d local_ij = this->grid->getLocalCoordinate();
    
    // initialize all data fields
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

    // begin solve
    double dt = this->getTimeStep();
    while (this->t <= this->t_end) {
        step_no++;
        this->t += dt;
        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        this->advanceSolution(dt,comm);
        
        MPI_Barrier(comm);
        if((step_no % write_every) == 0 || this->t >= this->t_end)
            this->write(comm);
    }
    
}

void FDHeatTransfer::advanceSolution(double dt,MPI_Comm comm)
{
    // share data
    for (auto p : this->data)
        FDModel::getDataFromNeighbors(getData(p.first), getSharedData(p.first), comm);
    
    // update boundary conditions
    this->updateBoundaryConditions(getSharedData("T"),comm);
    MPI_Barrier(comm);
    
    
    vector<double> &Up = getData("Tp");
    vector<double> &U = getData("T");
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
            {
                this->grid->stencilPoints(ij, getData(p.first), getSharedData(p.first), getStencil(p.first));
            }
            
            /// diffusion term
            Stencil T = getStencil("T");
            source = this->calculateDiffusion(T,getStencil("K"));
            /// advection term
            if(hasSource("u") || hasSource("v")|| hasSource("w"))
                source += this->calculateAdvection(T,getStencil("u"),getStencil("v"),getStencil("w"));
            
            if(hasSource("pressure"))
                source += this->calculateDeformationEnergy(T,getStencil("pressure"));
            
            idx = this->grid->idxFromCoord(i , j);

            Up[idx] =  U[idx] + dt * source;
            
        }
    }
//    size_t k = 0;
//    for(double ui : U)
//    {
//        printf("U[%d] = %f, Up[%k] = %f\n",k, ui,k, Up[k++]);
//    }
    U = Up;
}

void FDHeatTransfer::write(MPI_Comm comm)
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
    vector<double> &U = getData("T");
    
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

double FDHeatTransfer::getTimeStep()
{
    return MIN(this->hx*this->hx, this->hy*this->hy)*this->CFL;
}

