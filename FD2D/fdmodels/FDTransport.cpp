#include "FDTransport.h"
#include <cfloat>

FDTransport::FDTransport(string variable_name, int dim)
: FDModel(dim), primary_variable(variable_name)
{
  bcs = {nullptr,nullptr,nullptr,nullptr};
}

FDTransport::~FDTransport()
{
  
}

void FDTransport::initModel()
{
  // set local storage
  int local_els = this->grid->getNumberOfGridPoints();
  vector<double> zeros(local_els,0);
  
  this->data_manager->setData(primary_variable, zeros);
  this->data_manager->setData(primary_variable+"p", zeros);
}

void FDTransport::applyBoundaryConditions(MPI_Comm comm)
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

void FDTransport::updateBoundaryConditions(vector<vector<double>> &nbrs_data, MPI_Comm comm)
{
  if(!this->grid->isBoundaryGrid())
  {
    printf("NO BOUNDARY GRID!!!\n");
    return;
  }
  
  const vector<int> &nbrs = this->grid->getNeighbors();
  vector<double> &U = this->data_manager->getData(primary_variable);

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
  if (!this->data_manager->hasSource(initial_condition)) return;
  // want to get left and right values from neihbors
  int rank , nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  
  vector<double> &U = this->data_manager->getData(primary_variable);
  PhysicalSource *ic =  this->data_manager->getSource(initial_condition);
  double x,y;
  for (int idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++) {
    x = this->grid->getX(idx);
    y = this->grid->getY(idx);
    U[idx] = ic->operator()(x, y, 0, this->t_start);;
  }
  
  MPI_Barrier(comm);

}

double FDTransport::calculateAdvection(const FDUtils::Stencil &U,
                                       const FDUtils::Stencil &u,
                                       const FDUtils::Stencil &v,
                                       const FDUtils::Stencil &w)
{
  double advection = 0;
  
  if(this->dim >= 1)
  {
    double hx = this->grid->get_hx();
    if(u.O >= 0)
      advection += u.O * (U.O - U.W)/hx;
    else
      advection += u.O * (U.E - U.O)/hx;
  }
  
  if(this->dim >= 2)
  {
    double hy = this->grid->get_hy();
    if(v.O >= 0)
      advection += v.O * (U.O - U.S)/hy;
    else
      advection += v.O * (U.N - U.O)/hy;
  }
  
  if(this->dim >= 3)
  {
    double hz = this->grid->get_hz();
    if(w.O >= 0)
      advection += w.O * (U.O - U.B)/hz;
    else
      advection += w.O * (U.T - U.O)/hz;
  }
  
  advection += U.O * this->calculateDivergence(u,v,w);
  return -1.0 * advection;
}

void FDTransport::solve(MPI_Comm comm)
{
  int rank , nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  
  size_t step_no = 0;
  
  this->write(comm);
  
  // get new dt
  double dt;
  this->updateSources();
  while (this->t <= this->t_end) {
    
    // 0 update sources
    dt = this->getTimeStep(comm);
    
    // advance time step
    this->t += dt;
    step_no++;
    
    
    this->advanceSolution(dt, comm);
    
    // set solution
    this->data_manager->getData(primary_variable) = this->data_manager->getData(primary_variable+"p");
    MPI_Barrier(comm);
    if((step_no % write_every) == 0 || this->t >= this->t_end)
    {
      if(rank == 0) printf("completed step: %lu, time: %f[s]\n",step_no, this->t);
      this->write(comm);
    }
  }
  
}

void FDTransport::advanceSolution(double dt, MPI_Comm comm)
{
  // share data
  for (auto p : this->data_manager->availableData())
    FDModel::getDataFromNeighbors(this->data_manager->getData(p),
                                  this->data_manager->getSharedData(p), comm);
  
  // update boundary conditions
  this->updateBoundaryConditions(this->data_manager->getSharedData(primary_variable),comm);
  MPI_Barrier(comm);
  
  vector<double> &Up = this->data_manager->getData(primary_variable+"p");
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
        this->grid->stencilPoints(ij,
                                  this->data_manager->getData(p),
                                  this->data_manager->getSharedData(p),
                                  this->data_manager->getStencil(p));
      
      /// advection term
      if(this->data_manager->hasData("u") ||
         this->data_manager->hasData("v") ||
         this->data_manager->hasData("w"))
      {
        source = this->calculateAdvection(this->data_manager->getStencil(primary_variable),
                                          this->data_manager->getStencil("u"),
                                          this->data_manager->getStencil("v"),
                                          this->data_manager->getStencil("w"));
      }
      
      idx = this->grid->idxFromCoord(i , j);
      Up[idx] =  U[idx] + dt * source;
      
    }
  }
  //    U = Up;
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

double FDTransport::calculate_min_dt()
{
  vector<double> &u = this->data_manager->getData("u");
  vector<double> &v = this->data_manager->getData("v");
  vector<double> &w = this->data_manager->getData("w");
 
  double hx = this->grid->get_hx();
  double hy = this->grid->get_hy();
  double hz = this->grid->get_hz();
  double min_dt = DBL_MAX, dt = 0;
  double fact;
  for(size_t k = 0; k < this->grid->getNumberOfGridPoints(); k++)
  {
    fact =  fabs(u[k]/hx + v[k]/hy + w[k]/hz);
    dt = this->CFL /fact;
    if(dt < min_dt) min_dt = dt;
  }
  
  return min_dt;

}

double FDTransport::getTimeStep(MPI_Comm comm)
{
  int rank = MPI_Comm_rank(comm, &rank);
  double local_dt = calculate_min_dt();
  double dt = 1.0;
  MPI_Barrier(comm);
  MPI_Allreduce(&local_dt, &dt, 1, MPI_DOUBLE, MPI_MIN,comm);
  return fabs(dt) < 1e-8 ? 1e-8 : dt;
}

