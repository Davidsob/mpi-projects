#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <mpi.h>

#include <communicatorMap.h>
#include <LocalGrid.h>

#include <FDHeatTransfer.h>

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
    
    int divx = 15;
    int divy = 15;
    int Nx = divx+1;
    int Ny = divy+1;
    
    // Partition Global Grid
    FDHeatTransfer * model = new FDHeatTransfer(2);

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
                model->getGrid(), sizeof(LocalGrid),MPI_BYTE,
                0,MPI_COMM_WORLD);
    
    //  // scatter neighbor info to local grids
    vector<int> tmp(4,0);
    MPI_Scatter(all_neighbors.data(),4*sizeof(int),MPI_BYTE,
                tmp.data(), 4*sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);
    
    // set up model grid
    model->getGrid()->setNeighbors(tmp);
    MPI_Barrier(MPI_COMM_WORLD);

    // set up model
    // Courant number
    double CFL = 0.05;
    double t_end = 1.0;
    double hx = 1.0/(Nx-1.0) , hy = 1.0/(Ny - 1.0);
    
    model->setEndTime(t_end);
    model->setGridSpacing(hx,0);
    model->setGridSpacing(hy,1);
    model->setCFL(CFL);
    
    // set write data info
    model->setOuputFilePath(saveto + "/solutions/");
    model->setFileName("u.dat");
    model->setWriteEveryNthStep(20);
    MPI_Barrier(MPI_COMM_WORLD);
    
    // initialize solution vectors
    model->initModel();
    
    // set physical properties
    ConstantSource * unit  = new ConstantSource(1.0);
    model->setSource("rho", unit,true);
    model->setSource("K", unit,true);
    model->setSource("Cp", unit,true);
    
    // set boundary conditions
    double Tinit = 500;
    double T_amb = 300;
    double hc = -1.0;
    double insulation = 0;
    
    dirchletBoundaryCondition * fixed = new dirchletBoundaryCondition(Tinit);
    convectiveCooling * cooling_SN = new convectiveCooling(hc,T_amb,hy);
    convectiveCooling * cooling_E = new convectiveCooling(hc,T_amb,hx);
    model->setBoundaryConditions({fixed,cooling_SN,cooling_E,cooling_SN});
    
    // set velocity components
    ConstantSource * u = new ConstantSource(-10);
    ConstantSource * v = new ConstantSource(20);
    ConstantSource * w = new ConstantSource(0);
    model->setSource("u", u,true);
    model->setSource("v", v,true);
    model->setSource("w", w,true);
    
    // set initial condition
    ConstantSource * ic = new ConstantSource(Tinit);
    model->setSource("initial condition",ic,true);
    
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
    delete model;
    delete unit;
    
    delete u;
    delete v;
    delete w;
    
    delete fixed;
    delete cooling_E;
    delete cooling_SN;
    delete ic;
    
    MPI_Finalize();
}
