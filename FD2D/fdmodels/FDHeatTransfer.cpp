#include "FDHeatTransfer.h"

FDHeatTransfer::FDHeatTransfer(int dim)
: FDModel(dim)
{
    rho = nullptr;
    Cp = nullptr;
    K = nullptr;
    u = nullptr;
    v = nullptr;
    w = nullptr;
    p = nullptr;
    bcs = {nullptr,nullptr,nullptr,nullptr};
}

FDHeatTransfer::~FDHeatTransfer()
{
    //    if(rho) delete rho;
    //    if(Cp) delete Cp;
    //    if(K) delete Cp;
    //    if(u) delete u;
    //    if(v) delete v;
    //    if(w) delete w;
    //    if(p) delete p;
    //    for(boundaryCondition * bc : this->bcs)
    //      if(bc) delete bc;
}

void FDHeatTransfer::initModel()
{
    // set local storage
    int local_els = (grid->getRows()+1)*(grid->getCols()+1);
    vector<double> zeros(local_els,0);
    this->setU(zeros);
    this->setUp(zeros);
}

void FDHeatTransfer::applyBoundaryConditions(MPI_Comm comm)
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
                vector<double> boundary_u = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i), this->U);
                dynamic_cast<convectiveCooling *>(bcs[i])->updateBC(nbrs_data[i], boundary_u);
            }else{
                this->bcs[i]->updateBC(nbrs_data[i]);
            }
        }
    }
}

void FDHeatTransfer::applyInitialConditions(MPI_Comm comm)
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
    double dt = this->getTimeStep();
    
    this->write(comm);
    
    // get local informaion
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    Point2d local_ij = this->grid->getLocalCoordinate();
    
    // get
    double x,y;
    vector<double> Kij(this->U.size(),0);
    vector<double> u_ij(this->U.size(),0);
    vector<double> v_ij(this->U.size(),0);
    vector<double> w_ij(this->U.size(),0);
    vector<double> p_ij(this->U.size(),0);

    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            Kij[idx] = 1.0;
            if(this->u) u_ij[idx] = this->u->operator()(x, y, 0, dt);
            if(this->v) v_ij[idx] = this->v->operator()(x, y, 0, dt);
            if(this->w) w_ij[idx] = this->w->operator()(x, y, 0, dt);
            if(this->p) p_ij[idx] = this->p->operator()(x, y, 0, dt);
        }
    }
    
    MPI_Barrier(comm);
    
    // initialze stenice data
    Stencil T, Kmat, su, sv, sw, sp;

    vector<vector<double>> neighbors_data(4,{0});
    vector<vector<double>> neighbors_k(4,{0});
    vector<vector<double>> neighbors_u(4,{0});
    vector<vector<double>> neighbors_v(4,{0});
    vector<vector<double>> neighbors_w(4,{0});
    vector<vector<double>> neighbors_p(4,{0});
    Point2d ij{0,0};
    double source = 0;
    while (this->t <= this->t_end) {
        
        // advance time step
        this->t += dt;
        step_no++;
        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        // share data
        FDModel::getDataFromNeighbors(this->U,neighbors_data,comm);
        FDModel::getDataFromNeighbors(Kij,neighbors_k,comm);
        if(this->u) FDModel::getDataFromNeighbors(u_ij,neighbors_u,comm);
        if(this->v) FDModel::getDataFromNeighbors(v_ij,neighbors_v,comm);
        if(this->w) FDModel::getDataFromNeighbors(w_ij,neighbors_w,comm);
        if(this->p) FDModel::getDataFromNeighbors(p_ij,neighbors_p,comm);
        
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
                
                this->grid->stencilPoints(ij,this->U, neighbors_data, T);
                this->grid->stencilPoints(ij,Kij,neighbors_k, Kmat);
                if(this->u) this->grid->stencilPoints(ij, u_ij, neighbors_u, su);
                if(this->v) this->grid->stencilPoints(ij, v_ij, neighbors_v, sv);
                if(this->w) this->grid->stencilPoints(ij, w_ij, neighbors_w, sw);
                if(this->p) this->grid->stencilPoints(ij, p_ij, neighbors_p, sp);
                
                /// diffusion term
                source = this->calculateDiffusion(T,Kmat);
                /// advection term
                if(this->u || this->v || this->w)
                    source += this->calculateAdvection(T,su,sv,sw);
                
                if(this->p)
                    source += this->calculateDeformationEnergy(T,sp);
                
                this->Up[idx] =  this->U[idx] + dt * source;
                
            }
        }
        
        this->U = this->Up;
        MPI_Barrier(comm);
        this->write(comm);
    }
    
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

double FDHeatTransfer::getTimeStep()
{
    return MIN(this->hx*this->hx, this->hy*this->hy)*this->CFL;
}

