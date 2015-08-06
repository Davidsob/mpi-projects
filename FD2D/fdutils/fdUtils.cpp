#include "fdUtils.h"



namespace FDUtils {
  int idxFromCoord(int i, int j, int cols)
  { return i * cols + j;}
  
  // means
  double arithemeticMean(double x1, double x2){
    return (x1+x2)*0.5;
  }
  
  double harmonicMean(double x1, double x2){
    double f = (1.0/x1 + 1.0/x2)/2.0;
    return 1.0/f;
    
  }
  
  double geometricMean(double x1, double x2){
    return sqrt(x1*x2);
  }
}