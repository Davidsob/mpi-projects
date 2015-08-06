#pragma once

#include <random>
#include <stdint.h>
#include <vector>

class Random {

public:

  Random();
  ~Random();

  int                 getRandom(     int    min,          int    max         );
  std::vector<int>    getRandoms(    int    min,          int    max, int num);
  float               getRandom(     float  min,          float  max         );
  std::vector<float>  getRandoms(    float  min,          float  max, int num);
  double              getRandom(     double min,          double max         );
  std::vector<double> getRandoms(    double min,          double max, int num);

  float               getBreit(      float  mean = 0.0,   float gamma = 1.0  );
  float               getFDist(      float  m = 2.0,      float n = 2.0      );
  float               getExponential(float  lambda = 3.0                     );
  float               getGamma(      float  alpha = 2.0,  float beta = 2.0   );
  float               getWeiBull(    float  a = 2.0,      float b = 4.0      );
  float               getNormal(     float  mean = 0.0,   float sigma = 1.0  );
  float               getLogNormal(  float  mean = 0.0,   float sigma = 1.0  );
  float               getChiSq(      float  degrees = 3.0                    );
  float               getCauchy(     float  a = 5.0,      float b = 1.0      );
  float               getStudentT(   float  n = 10.0                         );
  bool                getBernoulli(  float  p = 0.5                          );
  int                 getPoisson(    float  mean = 4.0                       );
  int                 getBinomial(   int    t = 9,        float p = 0.5      );

private: // Methods

private: // Members

  std::random_device device;
  std::mt19937 engine;
};
