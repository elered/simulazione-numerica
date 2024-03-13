/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
using namespace std;

// This class contains functions for generating random numbers using the RANNYU algorithm
class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // Default constructor
  //Random();
  // Destructor
  ~Random();
  // Method to set the seed for the RNG
  void SetRandom(int * , int, int);
  // Method to save the seed to a file
  void SaveSeed();
  // Method to generate a random number in the range [0,1)
  Random();
  //Mio costruttore random
  double Rannyu(void);
  // Method to generate a random number in the range [min,max)
  double Rannyu(double min, double max);
  // Method to generate a random number with a Gaussian distribution
  double Gauss(double mean, double sigma);
  //Method to generate a random number with an exponential ditribution
  double Exp(double mean);
  //Method to generate a random number with a Cauchy-Lorentz distribution
  double CauchyLor (double gamma, double mu);

  double Cos();

  double Sin();
};

struct RandomWalkResult {
  
    vector<double> block_averages;
    vector<double> block_averages_squared;
    vector<double> block_errors;
};

RandomWalkResult random_walk_average(int num_steps, int num_trials, int num_blocks, double a, Random &rnd);

vector<double> mediarannyu(double numblocchi, double throws, Random rand);
vector<double> mediablocchi(double numblocchi, double throws, vector<double> media);
vector<double> mediablocchi2(double numblocchi, double throws, vector<double> dati);
vector<double> errore(double nblocchi, vector<double>mediablocchi, vector<double>mediablocchi2);
vector<double> integraleave(double xmin, double xmax, int throws, int blocchi, Random& rnd);
vector<double> integralehom(double xmin, double xmax, double fmax, Random& rnd, int throws, int blocchi);
void stampa(string nomeFile, double blocchi, const vector<double>integr, const vector<double>err_integr, const vector<double>integr_hom, const vector<double>err_hom);
//vector<double> random_walk_average(int num_steps, int num_trials, int num_blocks, double a, Random& rnd);
RandomWalkResult random_walk_average_cont(int num_steps, int num_trials, int num_blocks, double a, Random &rnd);

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
