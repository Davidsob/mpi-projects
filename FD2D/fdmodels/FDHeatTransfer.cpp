#include "FDHeatTransfer.h"

FDHeatTransfer::FDHeatTransfer(string variable_name, int dim)
: FDModel(dim), primary_variable(variable_name), coupled(false)
{
    bcs = {nullptr,nullptr,nullptr,nullptr};
}

FDHeatTransfer::~FDHeatTransfer()
{

}

void FDHeatTransfer::initModel()
{
    // set local storage
    int local_els = this->grid->getNumberOfGridPoints();
    vector<double> zeros(local_els,0);
    this->data_manager->setData(primary_variable, zeros);
    this->data_manager->setData(primary_variable+"p", zeros);
}

void FDHeatTransfer::applyBoundaryConditions(MPI_Comm comm)
{
    if(!this->grid->isBoundaryGrid())
    {
        printf("NO BOUNDARY GRID!!!\n");
        return;
    }
    
    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    vector<int> nbr = this->grid->getNeighbors();
    vector<double> &U = this->data_manager->getData(primary_variable);
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
    if(!this->grid->isBoundaryGrid())
    {
        printf("NO BOUNDARY GRID!!!\n");
        return;
    }
    
    const vector<int> &nbrs = this->grid->getNeighbors();
    
    for(size_t i = 0; i < nbrs.size(); i++)
    {
        if(nbrs[i] == -1)
        {
            if(this->bcs[i]->getName() == "convectiveCooling")
            {
                vector<double> boundary_u = this->grid->getBoundaryValues(FDUtils::GRID_DIRECTION(i),
                                                                          this->data_manager->getData(primary_variable));
                dynamic_cast<convectiveCooling *>(bcs[i])->updateBC(nbrs_data[i], boundary_u);
            }else{
                this->bcs[i]->updateBC(nbrs_data[i]);
            }
        }
    }
}

void FDHeatTransfer::applyInitialConditions(MPI_Comm comm)
{
    if (!this->data_manager->hasSource(initial_condition)) return;
    // want to get left and right values from neihbors
    int rank , nproc;
    MPI_Comm_size(comm, &nproc);
    MPI_Comm_rank(comm, &rank);
    
    // 1 set initial Temperature at 0
    vector<double> &U = this->data_manager->getData(primary_variable);
    PhysicalSource * ic = this->data_manager->getSource(initial_condition);
    double x,y;
    for (int idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++) {
        
        x = this->grid->getX(idx);
        y = this->grid->getY(idx);
        U[idx] = ic->operator()(x, y, 0, this->t_start);;
    }
    
    
    MPI_Barrier(comm);
}

double FDHeatTransfer::calculateDiffusion(const FDUtils::Stencil &T, const FDUtils::Stencil &K)
{
    double diff = 0.0;
   
    if(this->dim >= 1)
    {
        double hx = this->grid->get_hx();
        double rxsq = 1.0/pow(hx, 2);
        double d2x =  FDUtils::arithemeticMean(K.E,K.O)*(T.E - T.O) -
        FDUtils::arithemeticMean(K.W, K.O)*(T.O - T.W);
        diff += rxsq * d2x;
    }
    
    if(this->dim >= 2)
    {
        double hy = this->grid->get_hy();
        double rysq = 1.0/pow(hy, 2);
        double d2y = FDUtils::arithemeticMean(K.N,K.O)*(T.N - T.O) -
        FDUtils::arithemeticMean(K.S, K.O)*(T.O - T.S);
        diff += rysq * d2y;
    }
    
    if(this->dim == 3)
    {
        double hz = this->grid->get_hz();
        double rzsq = 1.0/pow(hz, 2);
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
        
        if(u.O >= 1)
            advection += u.O * (T.O - T.W)/this->grid->get_hx();
        else
            advection += u.O * (T.E - T.O)/this->grid->get_hx();
    }
    
    if(this->dim >= 2)
    {
        if(v.O >= 0)
            advection += v.O * (T.O - T.S)/this->grid->get_hy();
        else
            advection += v.O * (T.N - T.O)/this->grid->get_hy();
    }
    
    if(this->dim >= 3)
    {
        if(w.O >= 0)
            advection += w.O * (T.O - T.B)/this->grid->get_hz();
        else
            advection += w.O * (T.T - T.O)/this->grid->get_hz();
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
    
    size_t step_no = 0;
    // write initial solution
    this->write(comm);
    
    MPI_Barrier(comm);
    // begin solve
    double dt = this->getTimeStep(comm);
    if(rank == 0) printf("time step: %f[s]\n", dt);
    while (this->t <= this->t_end) {
        step_no++;
        this->updateSources();
        this->t += dt;
        if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
        
        this->advanceSolution(dt,comm);
      
        // set solution
        this->data_manager->getData(primary_variable) = this->data_manager->getData(primary_variable+"p");
      
        MPI_Barrier(comm);
        if((step_no % write_every) == 0 || this->t >= this->t_end)
            this->write(comm);
    }
    
}

void FDHeatTransfer::advanceSolution(double dt,MPI_Comm comm)
{
    
    // share data
    for (auto p : this->data_manager->availableData())
        FDModel::getDataFromNeighbors(this->data_manager->getData(p),
                                      this->data_manager->getSharedData(p), comm);
    
    // update boundary conditions
    this->updateBoundaryConditions(this->data_manager->getSharedData(primary_variable),comm);
    MPI_Barrier(comm);
    
    
    vector<double> &Up = this->data_manager->getData(primary_variable + "p");
    vector<double> &U = this->data_manager->getData(primary_variable);
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
    
            for (auto p : this->data_manager->availableData())
            {
                this->grid->stencilPoints(ij,
                                          this->data_manager->getData(p),
                                          this->data_manager->getSharedData(p),
                                          this->data_manager->getStencil(p));
            }
            
            /// diffusion term
            Stencil T = this->data_manager->getStencil(primary_variable);
            source = this->calculateDiffusion(T,this->data_manager->getStencil("K"));
            
            /// advection term
            if(this->data_manager->hasData("u") ||
               this->data_manager->hasData("v")||
               this->data_manager->hasData("w"))
            {
                source += this->calculateAdvection(
                                                   T,
                                                   this->data_manager->getStencil("u"),
                                                   this->data_manager->getStencil("v"),
                                                   this->data_manager->getStencil("w")
                                                   );
            }
            
            idx = this->grid->idxFromCoord(i , j);
            if(this->coupled)
            {
              if(this->data_manager->hasSource("pressure"))
                source += this->calculateDeformationEnergy(T,this->data_manager->getStencil("pressure"));
              
            }else{
              Up[idx] =  U[idx] + dt * source;
            }
            
        }
    }
  
//  if(!this->coupled) U = Up;
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
    
    vector<double> &U = this->data_manager->getData(primary_variable);
    
    for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++)
        file << this->grid->getX(idx) << "\t" << this->grid->getY(idx)  << "\t" << U[idx] <<"\n";
    
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

double FDHeatTransfer::getTimeStep(MPI_Comm comm)
{
    double hx = this->grid->get_hx();
    double hy = this->grid->get_hy();
    return MIN(hx*hx, hy*hy)*this->CFL;
}

