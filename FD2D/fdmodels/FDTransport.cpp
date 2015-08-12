#include "FDTransport.h"

FDTransport::FDTransport(int dim)
: FDModel(dim)
{
    u = nullptr;
    v = nullptr;
    w = nullptr;
    bcs = {nullptr,nullptr,nullptr,nullptr};
    maxU = 1;
}

FDTransport::~FDTransport()
{

    //    if(u) delete u;
    //    if(v) delete v;
    //    if(w) delete w;
    //    for(boundaryCondition * bc : this->bcs)
    //      if(bc) delete bc;
}

void FDTransport::initModel()
{
    // set local storage
    int local_els = (grid->getRows()+1)*(grid->getCols()+1);
    vector<double> zeros(local_els,0);
    this->setU(zeros);
    this->setUp(zeros);
}

void FDTransport::applyBoundaryConditions(MPI_Comm comm)
{
    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    vector<int> nbr = this->grid->getNeighbors();
    // set dirchlet
    if(nbr[0] == -1)
        for(i = 0; i <= local_rows; i++) this->bcs[0]->setBC(this->U[this->grid->idxFromCoord(i, 0)]);
    if(nbr[1] == -1)
        for(i = 0; i <= local_cols; i++) this->bcs[1]->setBC(this->U[this->grid->idxFromCoord(0, i)]);
    if(nbr[2] == -1)
        for(i = 0; i <= local_rows; i++) this->bcs[2]->setBC(this->U[this->grid->idxFromCoord(i, local_cols)]);
    if(nbr[3] == -1)
        for(i = 0; i <= local_cols; i++) this->bcs[3]->setBC(this->U[this->grid->idxFromCoord(local_rows, i)]);
    
    MPI_Barrier(comm);
}

void FDTransport::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
    if(!this->grid->isBoundaryGrid()) return;
    
    const vector<int> &nbrs = this->grid->getNeighbors();
    
    for(size_t i = 0; i < nbrs.size(); i++)
    {
        if(nbrs[i] == -1)
        {
            if(this->bcs[i]->getName() == "convectiveCooling")
            {
                vector<double> boundary_u = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i), this->U);
                dynamic_cast<convectiveCooling *>(bcs[i])->updateBC(nbrs_data[i], boundary_u);
            }else{
                this->bcs[i]->updateBC(nbrs_data[i]);
            }
        }
    }
}

void FDTransport::applyInitialConditions(MPI_Comm comm)
{
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
    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            this->U[idx] = this->initial_condition->operator()(x, y, 0, this->t_start);;
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
    vector<double> u_ij(this->U.size(),0);
    vector<double> v_ij(this->U.size(),0);
    vector<double> w_ij(this->U.size(),0);

    double localMaX = -1.0e24;
    double S;
    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            if(this->u) u_ij[idx] = this->u->operator()(x, y, 0, t);
            if(this->v) v_ij[idx] = this->v->operator()(x, y, 0, t);
            if(this->w) w_ij[idx] = this->w->operator()(x, y, 0, t);
            
            S = sqrt(pow(u_ij[idx], 2)+pow(v_ij[idx], 2)+pow(u_ij[idx], 3));
            if (S > localMaX) localMaX = S;
        }
    }
    
//    MPI_Barrier(comm);
//    printf("%d: %f\n", rank, localMaX);
//    MPI_Reduce(&localMaX,&maxU, sizeof(double)*1, MPI_DOUBLE,MPI_MAX,0,comm);
//    MPI_Bcast(&maxU, sizeof(double)*1,MPI_DOUBLE, 0, comm);
//    MPI_Barrier(comm);
//    printf("%d: %f\n", rank, maxU);
    // initialze stenice data
    Stencil rho, su, sv, sw, sp;

    vector<vector<double>> neighbors_data(4,{0});
    vector<vector<double>> neighbors_u(4,{0});
    vector<vector<double>> neighbors_v(4,{0});
    vector<vector<double>> neighbors_w(4,{0});
    Point2d ij{0,0};
    double source = 0;
    double dt = this->getTimeStep();
    printf("%d dt = %f\n",rank, dt);
    while (this->t <= this->t_end) {
        
        // advance time step
        this->t += dt;
        step_no++;
        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        // share data
        FDModel::getDataFromNeighbors(this->U,neighbors_data,comm);
        if(this->u) FDModel::getDataFromNeighbors(u_ij,neighbors_u,comm);
        if(this->v) FDModel::getDataFromNeighbors(v_ij,neighbors_v,comm);
        if(this->w) FDModel::getDataFromNeighbors(w_ij,neighbors_w,comm);
        
        // update boundary conditions
        this->updateBoundaryConditions(neighbors_data,comm);
        
        MPI_Barrier(comm);
        
        for(int i = 0; i <= local_rows; i++)
        {
            for(j = 0; j <= local_cols; j++)
            {
                ij.i = i;
                ij.j = j;
                idx = this->grid->idxFromCoord(i , j);
                
                this->grid->stencilPoints(ij,this->U, neighbors_data, rho);
                if(this->u) this->grid->stencilPoints(ij, u_ij, neighbors_u, su);
                if(this->v) this->grid->stencilPoints(ij, v_ij, neighbors_v, sv);
                if(this->w) this->grid->stencilPoints(ij, w_ij, neighbors_w, sw);
                
                /// advection term
                double source = 0;
                if(this->u || this->v || this->w)
                    source += this->calculateAdvection(rho,su,sv,sw);
                
                this->Up[idx] =  this->U[idx] + dt * source;
                
            }
        }
        
        this->U = this->Up;
        MPI_Barrier(comm);
        this->write(comm);
    }
    
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
    
    for (i = 0; i <= rows; i++) {
        for(j = 0; j <= cols; j++)
        {
            idx = this->grid->idxFromCoord(i, j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            file << x << "\t" << y << "\t" << this->U[idx] <<"\n";
        }
    }
    
    file.close();
    
    if(rank == 0) FDModel::file_count++;
}

double FDTransport::getTimeStep()
{
    return MIN(this->hx/maxU, this->hy/maxU)*this->CFL;
}

