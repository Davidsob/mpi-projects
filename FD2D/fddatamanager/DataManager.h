#pragma once

#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>

using namespace std;


class DataManager
{
public:
    virtual ~DataManager(){}
    
    virtual void setData(const string &name, vector<double> &d) = 0;
    virtual void setScalar(const string &name, double &s) = 0;
    
    virtual bool hasData(const string &name) = 0;
    virtual bool hasScalar(const string &name) = 0;
    
    virtual vector<double>& getData(const string & name) = 0;
    virtual double& getScalar(const string & name) = 0;
protected:
};


