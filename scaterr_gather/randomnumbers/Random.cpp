#include "Random.h"

#include <iostream>

Random::Random()
  : engine(device()){
}

Random::~Random() {

}

// ints
int Random::getRandom(int min, int max) {
  std::uniform_int_distribution<int> intDist(min, max);
  int number = intDist(engine);
  return number;
}

std::vector<int> Random::getRandoms(int min, int max, int num) {
  std::uniform_int_distribution<int> intDist(min, max);
  std::vector<int> numbers;
  for (int i = 0; i < num; ++i) {
    int number = intDist(engine);
    numbers.push_back(number);
  }
  return numbers;
}

// floats
float Random::getRandom(float min, float max) {
  std::uniform_real_distribution<float> floatDist(min, max);
  float number = floatDist(engine);
  return number;
}

std::vector<float> Random::getRandoms(float min, float max, int num) {
  std::uniform_real_distribution<float> floatDist(min, max);
  std::vector<float> numbers;
  for (int i = 0; i < num; ++i) {
    float number = floatDist(engine);
    numbers.push_back(number);
  }
  return numbers;
}

// doubles
double Random::getRandom(double min, double max) {
  std::uniform_real_distribution<double> doubleDist(min, max);
  double number = doubleDist(engine);
  return number;
}

std::vector<double> Random::getRandoms(double min, double max, int num) {
  std::uniform_real_distribution<double> doubleDist(min, max);
  std::vector<double> numbers;
  for (int i = 0; i < num; ++i) {
    double number = doubleDist(engine);
    numbers.push_back(number);
  }
  return numbers;
}

float Random::getBreit(float mean, float gamma) {
  float val, displ;
  float pi = std::atan(1)*4;
  val = 2*this->getRandom(0.0, 1.0) - 1;
  displ = 0.5*gamma*std::tan(val*(pi/2));
  return mean+displ;
}

float Random::getFDist(float m, float n) {
  std::fisher_f_distribution<float> dist(m,n);
  float val = dist(engine);
  return val;
}

bool Random::getBernoulli(float p) {
  // returns true or false based on probability p
  std::bernoulli_distribution dist(p);
  bool val = dist(engine);
  return val;
}

int Random::getPoisson(float mean) {
  std::poisson_distribution<int> dist(mean);
  int val = dist(engine);
  return val;
}

float Random::getGamma(float alpha, float beta) {
  std::gamma_distribution<float> dist(alpha, beta);
  float val = dist(engine);
  return val;
}

float Random::getWeiBull(float a, float b) {
  std::weibull_distribution<float> dist(a, b);
  float val = dist(engine);
  return val;
}

float Random::getNormal(float mean, float sigma) {
  std::normal_distribution<float> dist(mean, sigma);
  float val = dist(engine);
  return val;
}

float Random::getLogNormal(float mean, float sigma) {
  std::lognormal_distribution<float> dist(mean, sigma);
  float val = dist(engine);
  return val;
}

float Random::getChiSq(float degrees) {
  std::chi_squared_distribution<float> dist(degrees);
  float val = dist(engine);
  return val;
}

float Random::getCauchy(float a, float b) {
  std::cauchy_distribution<float> dist(a, b);
  float val = dist(engine);
  return val;
}

float Random::getStudentT(float n) {
  std::student_t_distribution<float> dist(n);
  float val = dist(engine);
  return val;
}

float Random::getExponential(float lambda) {
  std::exponential_distribution<float> dist(lambda);
  float val = dist(engine);
  return val;
}

int Random::getBinomial(int t, float p) {
  std::binomial_distribution<int> dist(t, p);
  int val = dist(engine);
  return val;
}
