//
//  FDModel.cpp
//  FD2D
//
//  Created by Bradley Davidson on 8/9/15.
//
//

#include "FDModel.h"

size_t FDModel::file_count = 0;

void FDModel::getDataFromNeighbors(vector<double> &data,
                                   vector<vector<double>> &nbrs_data,
                                   MPI_Comm comm)
{
    // pass receive depending on color
    
    if(this->grid->getColor() == 0)
    {
        // send
        this->grid->SendDataToNeighbor(WEST,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(SOUTH,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(EAST,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(NORTH,data,nbrs_data,comm);
        // receive
        this->grid->ReceiveDataFromNeighbor(WEST,nbrs_data[0],comm);
        this->grid->ReceiveDataFromNeighbor(SOUTH,nbrs_data[1],comm);
        this->grid->ReceiveDataFromNeighbor(EAST,nbrs_data[2],comm);
        this->grid->ReceiveDataFromNeighbor(NORTH,nbrs_data[3],comm);
        
    }else if(this->grid->getColor() == 1){
        
        // receive
        this->grid->ReceiveDataFromNeighbor(WEST,nbrs_data[0],comm);
        this->grid->ReceiveDataFromNeighbor(SOUTH,nbrs_data[1],comm);
        this->grid->ReceiveDataFromNeighbor(EAST,nbrs_data[2],comm);
        this->grid->ReceiveDataFromNeighbor(NORTH,nbrs_data[3],comm);
        
        // send
        this->grid->SendDataToNeighbor(WEST,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(SOUTH,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(EAST,data,nbrs_data,comm);
        this->grid->SendDataToNeighbor(NORTH,data,nbrs_data,comm);
    }
}

double FDModel::calculateDivergence(const FDUtils::Stencil &u,
                                    const FDUtils::Stencil &v,
                                    const FDUtils::Stencil &w)
{
    double div = 0;
    if(this->dim >= 1)
        div += (u.E - u.W)/this->grid->get_hx()/2.0;
    
    if(this->dim >= 2)
        div += (v.N - v.S)/this->grid->get_hy()/2.0;
    
    if(this->dim >= 3)
        div += (w.T - w.B)/this->grid->get_hz()/2.0;
    
    return div;
}

double FDModel::calculatePartialDerivative(const FDUtils::Stencil &U, int dim)
{
    if(this->dim == 0)
        return (U.E - U.W)/this->grid->get_hx()/2.0;
    
    if(this->dim == 1)
        return (U.N - U.S)/this->grid->get_hy()/2.0;
    
    if(this->dim == 2)
        return (U.T - U.B)/this->grid->get_hz()/2.0;
    
    return 0.0;
}

void FDModel::updateSources()
{
    // initialize all data fields
    double x, y;
    for (size_t idx = 0; idx < this->grid->getNumberOfGridPoints(); idx++) {
        
        x = this->grid->getX(idx);
        y = this->grid->getY(idx);
        
        for(auto p : this->data_manager->availableSources())
        {
            if(this->data_manager->hasData(p))
                this->data_manager->getData(p)[idx] =
                this->data_manager->getSource(p)->operator()(x, y, 0, this->t);
        }
    }
}