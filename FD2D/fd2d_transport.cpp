#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <mpi.h>

#include <communicatorMap.h>
#include <LocalGrid.h>

#include <FDTransport.h>
#include <FDDataManager.h>

#define MIN(A,B) A < B ? A : B


using namespace std;
using namespace FDSource;
using namespace FDGrid;
using namespace FDUtils;

static const string saveto = "/Users/LNLB/mpi-projects/FD2D";

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
    
    int divx = 40;
    int divy = 40;
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
    double CFL = 0.075;
    double t_end = 5.0;
    
    FDTransport * model = new FDTransport("density", 2);
    model->setGrid(grid);
    model->setDataManager(data);
    model->setEndTime(t_end);
    model->setCFL(CFL);
    
    // set write data info
    model->setOuputFilePath(saveto + "/solutions/");
    model->setFileName("u.dat");
    model->setWriteEveryNthStep(50);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // initialize solution vectors
    model->initModel();
    
    // set boundary conditions
    double insulation = 0;
    neumannBoundaryCondition * noflux_SN = new neumannBoundaryCondition(insulation, hy);
    neumannBoundaryCondition * noflux_WE = new neumannBoundaryCondition(insulation, hx);
    model->setBoundaryConditions({noflux_WE,noflux_SN,noflux_WE,noflux_SN});
    
    // set velocity components
    RadialSource2D * u = new RadialSource2D(0.5,0.5,4.0,0);
    RadialSource2D * v = new RadialSource2D(0.5,0.5,4.0,1);
    ConstantSource * w = new ConstantSource(0);
    data->setSource("u", u,true);
    data->setSource("v", v,true);
    data->setSource("w", w,true);
    
    // set initial condition
    GaussianSource * ic  = new GaussianSource(0.9,0.5,0.1,0.1,2.0);
    string ic_name = "IC";
    model->setICName(ic_name);
    data->setSource(ic_name, ic);
    
    // set IC
    if(rank == 0) printf("Applying initial conditions...\n");
    model->applyInitialConditions(MPI_COMM_WORLD);
    
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
    delete noflux_SN;
    delete noflux_WE;
    delete ic;
//    delete u;
//    delete v;
//    delete w;
    delete model;
    delete data;
    delete grid;
    MPI_Finalize();
}
