#pragma once
#include "DataManager.h"
#include <physicalSources.h>
#include <fdUtils.h>
#include <set>
using namespace FDSource;
using namespace FDUtils;

class FDDataManager : public DataManager
{
public:
    
    FDDataManager(size_t count = 0) : data_count(count){};
    
    void setDataCount(size_t count){ this->data_count = count;}
    virtual void setData(const string & name, vector<double> &d);
    virtual void setSource(const string & name, PhysicalSource * S, bool add_data = false);
    virtual void setScalar(const string & name, double &s);
    
    virtual bool hasData(const string & name);
    bool hasSource(const string & name);
    virtual bool hasScalar(const string & name);
    
    virtual vector<double>& getData(const string & name);
    vector<vector<double>>& getSharedData(const string & name);
    PhysicalSource * getSource(const string & name);
    FDUtils::Stencil & getStencil(const string & name);
    virtual double& getScalar(const string & name);

    set<const string>& availableData();
    set<const string>& availableSources();


private:

    void setSharedData(const string & name, vector<vector<double>> &sd);
    void setStencil(const string & name, FDUtils::Stencil &S);

    map<const string, vector<double>> data;
    map<string, vector<vector<double>> > shared_data; 
    map<const string, double> scalars;
    map<const string, PhysicalSource *> sources;
    map<const string, FDUtils::Stencil> stencils;
    
    size_t data_count;
    set<const string> registered_data, registered_sources;
    
};


