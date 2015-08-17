#include "FDCompressibleFlow.h"
#include "FDHeatTransfer.h"
#include "FDTransport.h"

FDCompressibleFlow::FDCompressibleFlow(vector<const string> variable_names, int dim)
: FDModel(dim), primary_variables(variable_names)
{
  FDHeatTransfer * ht = new FDHeatTransfer("temperature", dim);
  ht->setCoupled(true);
  FDTransport * mass = new FDTransport("density", dim);
  
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

void FDCompressibleFlow::setGrid(LocalGrid * _grid){
  this->grid = _grid;
  for(auto m : this->couplings) m.second->setGrid(_grid);
}

void FDCompressibleFlow::setDataManager(FDDataManager * dmngr){
  this->data_manager = dmngr;
  for(auto m : this->couplings) m.second->setDataManager(dmngr);
}

void FDCompressibleFlow::initModel()
{
  // set local storage
  int local_els = this->grid->getNumberOfGridPoints();
  vector<double> zeros(local_els,0);
  for(auto var : primary_variables)
  {
    this->data_manager->setData(var, zeros);
    this->data_manager->setData(var + "p", zeros);
  }
  
  this->data_manager->setData("divU", zeros);
  this->data_manager->setData("px", zeros);
  this->data_manager->setData("py", zeros);
  
  for(auto m : this->couplings) m.second->initModel();
}


void FDCompressibleFlow::applyBoundaryConditions(MPI_Comm comm)
{
  for(auto m : this->couplings)
  {
    printf("Applying bc for: %s\n", m.first.c_str());
    couplings[m.first]->applyBoundaryConditions(comm);
  }
  
  size_t i = 0;
  int local_rows = this->grid->getRows();
  int local_cols = this->grid->getCols();
  
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

void FDCompressibleFlow::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
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

void FDCompressibleFlow::applyInitialConditions(MPI_Comm comm)
{
  for(auto m : this->couplings) m.second->applyInitialConditions(comm);
  
  // want to get left and right values from neihbors
  int rank , nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  
  for(auto icpair : initial_conditions)
  {
    // 1 set initial Temperature at 0
    vector<double> &U = this->data_manager->getData(icpair.first);
    
    // get ic for var
    PhysicalSource * ic = this->data_manager->getSource(icpair.second);
    
    double x, y;
    for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++) {
      
      x = this->grid->getX(idx);
      y = this->grid->getY(idx);
      U[idx] = ic->operator()(x, y, 0, this->t_start);
    }
  }
  
  
  MPI_Barrier(comm);
}

double FDCompressibleFlow::calculateDiffusion(const FDUtils::Stencil &T, const FDUtils::Stencil &K)
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


double FDCompressibleFlow::calculateAdvection(const FDUtils::Stencil &T,
                                              const FDUtils::Stencil &u,
                                              const FDUtils::Stencil &v,
                                              const FDUtils::Stencil &w,
                                              const FDUtils::Stencil &div)
{
  double advection = 0;
  if(this->dim >= 1)
  {
    double hx = this->grid->get_hx();
    if(u.O >= 0)
      advection += u.O * (T.O - T.W)/hx;
    else
      advection += u.O * (T.E - T.O)/hx;
  }
  
  if(this->dim >= 2)
  {
    double hy = this->grid->get_hy();
    if(v.O >= 0)
      advection += v.O * (T.O - T.S)/hy;
    else
      advection += v.O * (T.N - T.O)/hy;
  }
  
  if(this->dim == 3)
  {
    double hz = this->grid->get_hz();
    if(w.O >= 0)
      advection += w.O * (T.O - T.B)/hz;
    else
      advection += w.O * (T.T - T.O)/hz;
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
  
  // 0. calculate pressure,
  //  momentum, and divU.
  this->updateSources();
  size_t step_no = 0;
  // write initial solution
  this->write(comm);
  
  MPI_Barrier(comm);
  
  // begin solve
  FDModel * mass = getModel("mass transport");
  FDModel * energy = getModel("energy");
  double dt = 0;
  while (this->t <= this->t_end) {
    
    dt = this->getTimeStep(comm);
    this->setTime(this->t + dt);
    this->updateSources();
    
    step_no++;
    if(rank == 0) printf("begin step: %lu, time: %f[s]\n",step_no, this->t);
    
    // 1. update mas transport
    mass->advanceSolution(dt, comm);
    
    this->updateData();
    // 2. update velocity field
    this->advanceSolution(dt,comm);
    
    // 3. update temperature
    energy->advanceSolution(dt, comm);
    
    // 4. set solutions
    this->data_manager->getData("density") = this->data_manager->getData("densityp");
    
    for(auto var : primary_variables)
      this->data_manager->getData(var) = this->data_manager->getData(var + "p");
    
    this->data_manager->getData("temperature") = this->data_manager->getData("temperaturep");
    
    // 5. write data to file
    MPI_Barrier(comm);
    if((step_no % write_every) == 0 || this->t >= this->t_end)
      this->write(comm);
    
  }
  
}

void FDCompressibleFlow::advanceSolution(double dt,MPI_Comm comm)
{
  
  
  // share data
  for (auto p : this->data_manager->availableData())
    FDModel::getDataFromNeighbors(this->data_manager->getData(p),
                                  this->data_manager->getSharedData(p), comm);
  
  // update boundary conditions
  this->updateBoundaryConditions(this->data_manager->getSharedData("T"),comm);
  MPI_Barrier(comm);
  
  
  vector<double> &Up = this->data_manager->getData("Up");
  vector<double> &U = this->data_manager->getData("U");
  
  vector<double> &Vp = this->data_manager->getData("Up");
  vector<double> &V = this->data_manager->getData("U");
  
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
      for (auto p : this->data_manager->availableData())
      {
        this->grid->stencilPoints(ij, this->data_manager->getData(p),
                                  this->data_manager->getSharedData(p),
                                  this->data_manager->getStencil(p));
      }
      
      // update u
      /// diffusion term
      Stencil u = this->data_manager->getStencil("u");
      Stencil v = this->data_manager->getStencil("v");
      Stencil w = this->data_manager->getStencil("w");
      Stencil P = this->data_manager->getStencil("pressure");
      Stencil px = this->data_manager->getStencil("px");
      Stencil py = this->data_manager->getStencil("py");
      Stencil divU = this->data_manager->getStencil("divU");
      Stencil mu = this->data_manager->getStencil("mu");
      
      
      /// diffusion term
      source = this->calculateDiffusion(u,mu);
      // volumetric expansion
      source += mu.O * this->calculatePartialDerivative(divU, 0)/3.0;
      // advection term
      source -= this->calculateAdvection(px,u,v,w,divU);
      // pressure term
      source -= this->calculatePartialDerivative(P, 0);
      // gravitational body force
      if(this->data_manager->hasSource("gx"))
        source += this->data_manager->getStencil("gx").O * this->data_manager->getStencil("density").O;
      
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
      if(this->data_manager->hasSource("gy"))
        source += this->data_manager->getStencil("gy").O * this->data_manager->getStencil("density").O;
      
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
  
  //    printf("Data to write:\n");
  //    for(auto s : this->data_manager->availableData())
  //        printf("data soure: %s\n", s.c_str());
  
  vector<double> &U = this->data_manager->getData(primary_variables[0]);
  vector<double> &V = this->data_manager->getData(primary_variables[1]);
  vector<double> &W = this->data_manager->getData("w");
  vector<double> &P = this->data_manager->getData("pressure");
  vector<double> &T = this->data_manager->getData("temperature");
  vector<double> &rho = this->data_manager->getData("density");
  vector<double> &mu = this->data_manager->getData("mu");
  
  
  double x,y;
  for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++)
  {
    x = this->grid->getX(idx);
    y = this->grid->getY(idx);
    file << x << "\t" << y;
    file << "\t" << U[idx] <<"\t" << V[idx] << "\t" << W[idx] << "\t" << P[idx];
    file << "\t" << T[idx] <<"\t" << rho[idx] << "\t" << mu[idx] << "\n";
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
  vector<double> &rho = this->data_manager->getData("density");
  vector<double> &T = this->data_manager->getData("temperature");
  vector<double> &p = this->data_manager->getData("pressure");
  
  for(size_t idx = 0; idx < rho.size(); idx++){
    p[idx] = rho[idx] * R * T[idx];
  }
  
}

void FDCompressibleFlow::calculateGridDivergence()
{
  vector<double> &div = this->data_manager->getData("divU");
  size_t idx;
  Point2d ij;
  for (int i = 0; i <= this->grid->getRows(); i++) {
    for(int j = 0; j <= this->grid->getCols(); j++)
    {
      ij.i = i; ij.j = j;
      idx = this->grid->idxFromCoord(i , j);
      
      for(auto s : primary_variables)
        this->grid->stencilPoints(ij,
                                  this->data_manager->getData(s),
                                  this->data_manager->getSharedData(s),
                                  this->data_manager->getStencil(s));
      
      div[idx] = this->calculateDivergence(this->data_manager->getStencil("u"),
                                           this->data_manager->getStencil("v"),
                                           this->data_manager->getStencil("w"));
    }
  }
}

void FDCompressibleFlow::calculateMomentum()
{
  vector<double> &rho = this->data_manager->getData("density");
  vector<double> &U = this->data_manager->getData("u");
  vector<double> &V = this->data_manager->getData("v");
  vector<double> &px = this->data_manager->getData("px");
  vector<double> &py = this->data_manager->getData("py");
  
  size_t k = 0;
  for(const double &rho_i : rho){
    px[k] = rho_i * U[k];
    py[k] = rho_i * V[k];
    k++;
  }
}

double FDCompressibleFlow::getTimeStep(MPI_Comm comm)
{
  FDModel * m = getModel("mass transport");
  return dynamic_cast<FDTransport *>(m)->getTimeStep(comm);
}

