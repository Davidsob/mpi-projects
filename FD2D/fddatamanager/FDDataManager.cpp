#include "FDDataManager.h"


void FDDataManager::setData(const string & name, vector<double> &d)
{
    registered_data.insert(name);
    this->data.insert(pair<const string, vector<double>>(name, d));
    
    // add shared data by same name
    vector<vector<double>> tmp(4, {0});
    this->setSharedData(name, tmp);
    
    // add stencil by same name
    Stencil S;
    this->setStencil(name, S);
}

void FDDataManager::setSource(const string & name,PhysicalSource * S,  bool add_data)
{
    registered_sources.insert(name);
    this->sources.insert(pair<const string, PhysicalSource *>(name, S));
    if(add_data)
    {
        // depending on use may or may not want
        // an associated data cache
        vector<double> tmp(this->data_count,0);
        this->setData(name, tmp);
    }
}

void FDDataManager::setSharedData(const string & name, vector<vector<double>> &sd)
{
    this->shared_data.insert(pair<const string, vector<vector<double>> >(name, sd));
}
void FDDataManager::setStencil(const string & name, FDUtils::Stencil &S)
{
    this->stencils.insert(pair<const string, FDUtils::Stencil>(name, S));
}

void FDDataManager::setScalar(const string & name, double &s)
{
    this->scalars.insert(pair<const string, double>(name, s));
}

bool FDDataManager::hasData(const string & name)
{
    return this->data.find(name) != this->data.end();
}
bool FDDataManager::hasSource(const string & name)
{
    return this->sources.find(name) != this->sources.end();
}
bool FDDataManager::hasScalar(const string & name)
{
    return this->scalars.find(name) != this->scalars.end();
}

vector<double>& FDDataManager::getData(const string & name)
{
    return this->data[name];
}

vector<vector<double>>& FDDataManager::getSharedData(const string & name)
{
    return this->shared_data[name];
}

PhysicalSource * FDDataManager::getSource(const string & name)
{
    return this->sources[name];
}

FDUtils::Stencil & FDDataManager::getStencil(const string & name)
{
    return this->stencils[name];
}

double& FDDataManager::getScalar(const string & name)
{
    return this->scalars[name];
}

set<const string>& FDDataManager::availableData()
{
    return registered_data;
}

set<const string>& FDDataManager::availableSources()
{
    return registered_sources;
}




