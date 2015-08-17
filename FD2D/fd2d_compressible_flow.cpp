#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <mpi.h>

#include <communicatorMap.h>
#include <LocalGrid.h>

#include <FDCompressibleFlow.h>
#include <FDTransport.h>
#include <FDHeatTransfer.h>

#define MIN(A,B) A < B ? A : B


using namespace std;
using namespace FDSource;
using namespace FDGrid;
using namespace FDUtils;

static const string saveto = "/Users/davidson/mpi-projects/FD2D";

void writeGridPartition(const vector<LocalGrid> & parts)
{
    ofstream file;
    file.open(saveto + string("/grids/graph.dat"));
    for(const LocalGrid & p : parts) p.write(file);
    file.close();
}

void getWorldInfo(MPI_Comm comm, int &r, int &ws)
{
    MPI_Comm_rank(comm,&r);
    MPI_Comm_size(comm,&ws);
}


int main(int argc, char ** argv)
{
    MPI_Init(NULL,NULL);
    
    // start process timer
    double timer = 0;
    timer -= MPI_Wtime();
    
    int rank, nproc;
    getWorldInfo(MPI_COMM_WORLD,rank,nproc);
    
    int divx = 20;
    int divy = 20;
    int Nx = divx+1;
    int Ny = divy+1;
    
    // Partition Global Grid
    LocalGrid * grid = new LocalGrid;

    vector<LocalGrid> mesh;
    vector<int> all_neighbors;
    vector<LocalGrid>grids;
    if(rank == 0)
    {
        CommunicatorMap map(nproc);
        
        for(int i = 0; i < nproc; i++)
        {
            LocalGrid lg;
            map.setLocalGrid(Nx, Ny, i, &lg);
            grids.push_back(lg);
        }
        
        // scatter nbrs
        all_neighbors = map.stackNbrs();
        
        // print the grid partition
        writeGridPartition(grids);
    }
    
    //  // scatter local grids
    MPI_Scatter(grids.data(),sizeof(LocalGrid),MPI_BYTE,
                grid, sizeof(LocalGrid),MPI_BYTE,
                0,MPI_COMM_WORLD);
    
    //  // scatter neighbor info to local grids
    vector<int> tmp(4,0);
    MPI_Scatter(all_neighbors.data(),4*sizeof(int),MPI_BYTE,
                tmp.data(), 4*sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);
    
    // finish set up of grid
    double hx = 1.0/(Nx - 1.0) , hy = 1.0/(Ny - 1.0);
    grid->setCellIncrements(hx, hy); // set spacing
    grid->initCoordinates(); // set coordinates
    grid->setNeighbors(tmp); // set neihbors
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // set up data manager
    FDDataManager * data = new FDDataManager(grid->getNumberOfGridPoints());

    // set up model
    // Courant number
    double CFL = 0.05;
    double t_end = 1.0;
    
    FDCompressibleFlow * model = new FDCompressibleFlow({"u", "v"}, 2);
    model->setGrid(grid);
    model->setDataManager(data);
    model->setEndTime(t_end);
    model->setCFL(CFL);
    
    // set write data info
    model->setOuputFilePath(saveto + "/solutions/");
    model->setFileName("u.dat");
    model->setWriteEveryNthStep(20);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // initialize solution vectors
    model->initModel();
    
    // set material properties
    ConstantSource * unit  = new ConstantSource(1.0);
    ConstantSource * zero = new ConstantSource(0.0);
    data->setSource("mu", unit,true);
    data->setSource("pressure", zero,true);
    
    // set boundary conditions for flow model
    dirchletBoundaryCondition * u_in = new dirchletBoundaryCondition(1);
    dirchletBoundaryCondition * wall = new dirchletBoundaryCondition(0);
    neumannBoundaryCondition * outflow_x = new neumannBoundaryCondition(0,hx);
    neumannBoundaryCondition * outflow_y = new neumannBoundaryCondition(0,hy);
    
    model->setBoundaryConditions("u", {u_in,wall,outflow_x,wall});
    model->setBoundaryConditions("v", {wall,wall,outflow_x,wall});
    
    // set velocity components (as a source)
//    RadialSource2D * u = new RadialSource2D(0.5, 0.5, 4.0, 0);
//    RadialSource2D * v = new RadialSource2D(0.5, 0.5, 4.0, 1);
    ConstantSource * w = new ConstantSource(0);
//    data->setSource("u", u);
//    data->setSource("v", v);
    data->setSource("w", w, true);
  
    // set initial condition on velocity
    ConstantSource * ic = new ConstantSource(0);
    model->addICName("u", "ic u");
    model->addICName("v", "ic v");
    data->setSource("ic u",ic);
    data->setSource("ic v",ic);
    
    // set up mass transport model
    FDTransport * mass = dynamic_cast<FDTransport *>(model->getModel("mass transport"));
    mass->setBoundaryConditions({u_in, outflow_y, outflow_x, outflow_y});
    
    // set up intial conditions
    GaussianSource * gauss = new GaussianSource(0.6,0.6,0.1,0.1,5.0);
    string mass_ic_name = "mass IC";
    mass->setICName(mass_ic_name);
    data->setSource(mass_ic_name, unit);
    
    // set up energy model
    FDHeatTransfer * heat = dynamic_cast<FDHeatTransfer *>(model->getModel("energy"));
    ConstantSource * thermal_conduct  = new ConstantSource(100.0);
    data->setSource("K", thermal_conduct,true);
    heat->setBoundaryConditions({outflow_x, outflow_y, outflow_x, outflow_y});
    
    // set initial condition
    ConstantSource * Tinit = new ConstantSource(300);
    string heat_ic_name = "heat IC";
    heat->setICName(heat_ic_name);
    data->setSource(heat_ic_name, Tinit);
    
    // set IC
    if(rank == 0) printf("Applying initial conditions...\n");
    model->applyInitialConditions(MPI_COMM_WORLD);
    //    setIC(u,hx,hy,C,dynamic_cast<PhysicalSource *>(H), my_grid, MPI_COMM_WORLD);
    
    // set BC
    if(rank == 0) printf("Applying boundary conditions...\n");
    model->applyBoundaryConditions(MPI_COMM_WORLD);
    
    // Solve
    if(rank == 0) printf("Solving...\n");
    model->solve(MPI_COMM_WORLD);
    
    // write file count
    if(rank == 0)
    {
        ofstream file;
        file.open(saveto + "/solutions/model_info.dat");
        file << FDModel::getFileCount() << "\n";
        file << nproc << "\n";
        file << model->getCFL() << "\n";
        file.close();
        
        timer += MPI_Wtime();
        printf("Simulation complete: time/proc = %f[s]\n",timer/nproc);
        printf("Simulation complete....\n");
    }
    
    // clean up
  delete thermal_conduct;
    delete unit;
    delete zero;
//    delete u;
//    delete v;
    delete w;
    delete u_in;
    delete wall;
    delete outflow_x;
    delete outflow_y;
    delete ic;
    delete gauss;
    delete Tinit;
    delete model;
    delete data;
    delete grid;
    
    MPI_Finalize();
}
