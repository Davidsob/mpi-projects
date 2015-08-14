#include "FDCompressibleFlow.h"
#include "FDHeatTransfer.h"
#include "FDTransport.h"

FDCompressibleFlow::FDCompressibleFlow(int dim)
: FDModel(dim)
{
    FDHeatTransfer * ht = new FDHeatTransfer(dim);
    FDTransport * mass = new FDTransport(dim);
    
    this->addModel("energy", ht);
    this->addModel("mass transport", mass);
}

FDCompressibleFlow::~FDCompressibleFlow()
{
    FDModel * heat = getModel("energy");
    delete dynamic_cast<FDHeatTransfer *>(heat);
    
    FDModel * mass = getModel("mass transport");
    delete dynamic_cast<FDTransport *>(mass);
}

void FDCompressibleFlow::initModel()
{
    // set local storage
    int local_els = (grid->getRows()+1)*(grid->getCols()+1);
    vector<double> zeros(local_els,0);
    this->setData("U", zeros);
    this->setData("Up", zeros);
    
    this->setData("V", zeros);
    this->setData("Vp", zeros);
    
    this->setData("W", zeros);
    this->setData("Wp", zeros);
    
    this->setData("divU", zeros);
    this->setData("px", zeros);
    this->setData("py", zeros);
}

void FDCompressibleFlow::applyBoundaryConditions(MPI_Comm comm)
{
    size_t i = 0;
    int local_rows = this->grid->getRows();
    int local_cols = this->grid->getCols();
    
    vector<int> nbr = this->grid->getNeighbors();
    
    // set dirchlet
    for(auto bc : this->bcs)
    {
        vector<double> &Ui = this->getData(bc.first);
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

void FDCompressibleFlow::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
    if(!this->grid->isBoundaryGrid()) return;
    
    const vector<int> &nbrs = this->grid->getNeighbors();
    
    for(auto bc : this->bcs)
    {
        vector<vector<double>> & nbr_data =  getSharedData(bc.first);
        vector<double> &Ui = getData(bc.first);
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

void FDCompressibleFlow::applyInitialConditions(MPI_Comm comm)
{
    this->getModel("mass transfer")->applyInitialConditions(comm);
    this->getModel("energy")->applyInitialConditions(comm);
    
    if (!hasSource("initial condition u") || !hasSource("initial condition v")) return;
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
    vector<double> &U = this->getData("U");
    vector<double> &V = this->getData("V");
    
    PhysicalSource * ic_u = getSource("initial condition u");
    PhysicalSource * ic_v = getSource("initial condition v");
    for (i = 0; i <= local_rows; i++) {
        for(j = 0; j <= local_cols; j++)
        {
            idx = this->grid->idxFromCoord(i,j);
            x = (local_ij.j + j)*this->hx;
            y = (local_ij.i + i)*this->hy;
            U[idx] = ic_u->operator()(x, y, 0, this->t_start);
            V[idx] = ic_v->operator()(x, y, 0, this->t_start);
        }
    }
    
    MPI_Barrier(comm);
}

double FDCompressibleFlow::calculateDiffusion(const FDUtils::Stencil &T, const FDUtils::Stencil &K)
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


double FDCompressibleFlow::calculateAdvection(const FDUtils::Stencil &T,
                                              const FDUtils::Stencil &u,
                                              const FDUtils::Stencil &v,
                                              const FDUtils::Stencil &w,
                                              const FDUtils::Stencil &div)
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
    
    advection += T.O * div.O;
    return advection;
}


double FDCompressibleFlow::calculateDeformationEnergy(const FDUtils::Stencil &p,const FDUtils::Stencil &u)
{
    double deformation = 0;
    return deformation;
}

void FDCompressibleFlow::solve(MPI_Comm comm)
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
        
        // 1. update mas transport
        getModel("mass transport")->advanceSolution(dt, comm);
        
        // 2. update velocity field
        this->advanceSolution(dt,comm);
        
        // 3. update temperature
        getModel("energy")->advanceSolution(dt, comm);
        
        MPI_Barrier(comm);
        if((step_no % write_every) == 0 || this->t >= this->t_end)
            this->write(comm);
    }
    
}

void FDCompressibleFlow::advanceSolution(double dt,MPI_Comm comm)
{
    // 0. calculate pressure,
    //  momentum, and divU.
    this->updateData();
    
    // share data
    for (auto p : this->data)
        FDModel::getDataFromNeighbors(getData(p.first), getSharedData(p.first), comm);
    
    // update boundary conditions
    this->updateBoundaryConditions(getSharedData("T"),comm);
    MPI_Barrier(comm);
    
    
    vector<double> &Up = getData("Up");
    vector<double> &U = getData("U");
    
    vector<double> &Vp = getData("Up");
    vector<double> &V = getData("U");
    
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
            idx = this->grid->idxFromCoord(i , j);

            // set all stencils
            for (auto p : this->stencils)
            {
                this->grid->stencilPoints(ij, getData(p.first), getSharedData(p.first), getStencil(p.first));
            }
            
            // update u
            /// diffusion term
            Stencil u = getStencil("U");
            Stencil v = getStencil("V");
            Stencil w = getStencil("W");
            Stencil P = getStencil("pressure");
            Stencil px = getStencil("px");
            Stencil py = getStencil("py");
            Stencil divU = getStencil("divU");
            Stencil mu = getStencil("mu");
            
            
            /// diffusion term
            source = this->calculateDiffusion(u,mu);
            // volumetric expansion
            source += mu.O * this->calculatePartialDerivative(divU, 0)/3.0;
            // advection term
            source -= this->calculateAdvection(px,u,v,w,divU);
            // pressure term
            source -= this->calculatePartialDerivative(P, 0);
             // gravitational body force
            if(this->hasSource("gx"))
                source += getStencil("gx").O * getStencil("density").O;
            
            Up[idx] =  U[idx] + dt * source;
            
            
            // update v
            /// diffusion term
            source = this->calculateDiffusion(v,mu);
            // volumetric expansion
            source += mu.O * this->calculatePartialDerivative(divU, 1)/3.0;
            // advection term
            source -= this->calculateAdvection(py,u,v,w,divU);
            // pressure term
            source -= this->calculatePartialDerivative(P, 1);
            // gravitational body force
            if(this->hasSource("gy"))
                source += getStencil("gy").O * getStencil("density").O;
            
            Vp[idx] =  V[idx] + dt * source;
            
        }
    }

    U = Up;
    V = Vp;
}

void FDCompressibleFlow::write(MPI_Comm comm)
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

void FDCompressibleFlow::updateData()
{
    this->calculatePressure();
    this->calculateGridDivergence();
    this->calculateMomentum();
}

void FDCompressibleFlow::calculatePressure()
{
    static const double R = 461.5;
    vector<double> &rho = this->getData("density");
    vector<double> &T = this->getData("temperature");
    vector<double> &p = this->getData("pressure");
    size_t k = 0;
    for(double &pi : p){
        pi = rho[k] * R * T[k];
        k++;
    }
    
}

void FDCompressibleFlow::calculateGridDivergence()
{
    vector<double> &div = this->getData("divU");
    size_t idx;
    Point2d ij;
    vector<string>vars = {"U", "V", "W"};
    for (int i = 0; i <= this->grid->getRows(); i++) {
        for(int j = 0; j <= this->grid->getCols(); j++)
        {
            ij.i = i; ij.j = j;
            idx = this->grid->idxFromCoord(i , j);
            
            for(const string &s : vars)
                this->grid->stencilPoints(ij, getData(s), getSharedData(s), getStencil(s));
            
            div[idx] = this->calculateDivergence(getStencil("U"), getStencil("V"), getStencil("W"));
        }
    }
}

void FDCompressibleFlow::calculateMomentum()
{
    vector<double> &rho = this->getData("density");
    vector<double> &U = this->getData("U");
    vector<double> &V = this->getData("V");
    vector<double> &px = this->getData("px");
    vector<double> &py = this->getData("py");
    
    size_t k = 0;
    for(const double &rho_i : rho){
        px[k] = rho_i * U[k];
        py[k] = rho_i * V[k];
        k++;
    }
}

double FDCompressibleFlow::getTimeStep()
{
    return MIN(this->hx*this->hx, this->hy*this->hy)*this->CFL;
}

