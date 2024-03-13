/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "random.h"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

// Random :: Random(){}
//  Default constructor, does not perform any action

Random ::~Random() {}
// Default destructor, does not perform any action

void Random ::SaveSeed() {
  // This function saves the current state of the random number generator to a
  // file "seed.out"
  ofstream WriteSeed;
  WriteSeed.open("seed.out");
  if (WriteSeed.is_open()) {
    WriteSeed << "RANDOMSEED	" << l1 << " " << l2 << " " << l3 << " " << l4
              << endl;
    ;
  } else
    cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

// mio costruttore random
Random ::Random() {
  int seed[4];
  int p1, p2;
  ifstream Primes("Primes");
  if (Primes.is_open()) {
    Primes >> p1 >> p2;
  } else
    cerr << "PROBLEM: Unable to open Primes" << endl;
  Primes.close();

  ifstream input("seed.in");
  string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        SetRandom(seed, p1, p2);
      }
    }
    input.close();
  } else
    cerr << "PROBLEM: Unable to open seed.in" << endl;
  SaveSeed();
}

double Random ::Gauss(double mean, double sigma) {
  // This function generates a random number from a Gaussian distribution with
  // given mean and sigma
  double s = Rannyu();
  double t = Rannyu();
  double x = sqrt(-2. * log(1. - s)) * cos(2. * M_PI * t);
  return mean + x * sigma;
}

double Random ::Rannyu(double min, double max) {
  // This function generates a random number in the range [min, max)
  return min + (max - min) * Rannyu();
}

double Random ::Rannyu(void) {
  // This function generates a random number in the range [0,1)
  const double twom12 = 0.000244140625;
  int i1, i2, i3, i4;
  double r;

  i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
  i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
  i3 = l3 * m4 + l4 * m3 + n3;
  i4 = l4 * m4 + n4;
  l4 = i4 % 4096;
  i3 = i3 + i4 / 4096;
  l3 = i3 % 4096;
  i2 = i2 + i3 / 4096;
  l2 = i2 % 4096;
  l1 = (i1 + i2 / 4096) % 4096;
  r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * (l4))));

  return r;
}

double Random ::Exp(double mean) {
  double y = Rannyu();
  return -(1 / mean) * log(1. - y);
}

double Random ::Cos() {
  double x = Rannyu(0, pow(10, 5));
  return cos(x);
}

double Random ::Sin() {
  double x = Rannyu(0, M_PI);
  return sin(x);
}

double Random ::CauchyLor(double gamma, double mu) {
  double x = Rannyu();

  return mu + gamma * tan(M_PI * (x - 0.5));
}

void Random ::SetRandom(int *s, int p1, int p2) {
  // This function sets the seed and parameters of the random number generator
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

vector<double> mediarannyu(double numblocchi, double throws, Random rand) {
  double L = throws / numblocchi;  // Number of throws in each block, please use
                                   // for M a multiple of N
  double sum = 0;                  // somme parziali
  vector<double> vec_media(numblocchi);
  vector<double> sumprog(numblocchi);

  for (int k = 0; k < numblocchi; k++) {
    sum = 0;

    for (int i = 0; i < L; i++) {
      sum += rand.Rannyu();
    }
    vec_media[k] = sum / L;
  }

  return vec_media;
}

vector<double> mediablocchi(double numblocchi, double throws, const vector<double> media) {
  vector<double> sumprog(numblocchi);

  for (int i = 0; i < numblocchi; i++) {
    for (int j = 0; j <= i; j++) {
      sumprog[i] += media[j];
    }
    sumprog[i] /= (i + 1);
  }

  return sumprog;
}

vector<double> mediablocchi2(double numblocchi, double throws, const vector<double> dati) {
  vector<double> vec_media2(numblocchi);
  vector<double> sumprog2(numblocchi);

  for (int k = 0; k < numblocchi; k++) {
    vec_media2[k] = dati[k] * dati[k];
  }

  for (int i = 0; i < numblocchi; i++) {
    for (int j = 0; j <= i; j++) {
      sumprog2[i] += vec_media2[j];
    }
    sumprog2[i] /= (i + 1);
  }

  return sumprog2;
}

vector<double> errore(double nblocchi, vector<double> mediablocchi, vector<double> mediablocchi2) {
  vector<double> var(nblocchi);
  for (int i = 0; i < nblocchi; i++) {
    if (i == 0) {
      var[i] = 0;
    } else {
      var[i] = sqrt((mediablocchi2[i] - pow(mediablocchi[i], 2)) / (i));
    }
  }

  return var;
}

vector<double> integraleave(double xmin, double xmax, int throws, int blocchi,Random &rnd) {
  double sum = 0;
  double L = throws / blocchi;
  double x = 0;
  vector<double> media_int(blocchi);

  for (int k = 0; k < blocchi; k++) {
    sum = 0;

    for (int i = 0; i < L; i++) {
      x = rnd.Rannyu(xmin, xmax);
      sum += (M_PI / 2) * cos((M_PI * x) / 2);
    }
    media_int[k] = (xmax - xmin) * (sum / L);
  }

  return media_int;
};

vector<double> integralehom(double xmin, double xmax, double fmax, Random &rnd,int throws, int blocchi) {
  double x = 0;

  double y = 0;

  double gx = 0;

  double L = throws / blocchi;

  vector<double> media_int(blocchi);

  for (int i = 0; i < blocchi; i++) {
    gx = 0;

    for (int j = 0; j < L; j++) {
      x = rnd.Rannyu(xmin, xmax);

      y = rnd.Rannyu(0, fmax);

      if (y < ((3 / M_PI) * (1 - (x * x)))) {
        gx += (M_PI / 3) * ((cos((M_PI * x) / 2)) / (1 - (x * x)));

      } else {
        j--;
      }
    }

    media_int[i] = gx / L;
  }

  return media_int;
};

void stampa(string nomeFile, double blocchi, const vector<double> integr,const vector<double> err_integr, const vector<double> integr_hom,const vector<double> err_hom) {
  fstream foutput;                   // foutput è il nome che do alla variabile
  foutput.open(nomeFile, ios::out);  // risultato.dat è il nome del file
  if (foutput.good())  // controllo che funzioni e gli faccio stampare,
                       // altrimenti gli faccio dare errore
  {
    for (double i = 0; i < blocchi; i++) {
      foutput << i + 1 << " " << integr[i] << " " << err_integr[i] << " "
              << integr_hom[i] << " " << err_hom[i] << endl;
    }
  } else
    cout << "Errore: file is not good!" << endl;

  foutput.close();  // chiudo il file
}


RandomWalkResult random_walk_average(int num_steps, int num_trials, int num_blocks, double a, Random &rnd) {
  vector<double> block_averages(num_steps, 0.0);
  vector<double> block_averages_squared(num_steps, 0.0);
  vector<double> block_errors(num_steps, 0.0);
  enum Direction { X, Y, Z };
  double sign = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double distance= 0;
  double block_sum = 0;

  // Itera su ogni passo del random walk
  for (int step = 1; step <= num_steps; step++) {
    // Itera su ogni blocco
    for (int block = 0; block < num_blocks; block++) {
      block_sum = 0;

      // Itera su ogni lancio nel blocco
      for (int trial = 0; trial < num_trials / num_blocks; trial++) {
        x = 0;
        y = 0;
        z = 0;

        // Esegue il random walk 3D
        for (int s = 0; s < step; s++) {
          Direction direction = static_cast<Direction>(rnd.Rannyu(0, 3));
          sign = rnd.Rannyu(0, 1);

          if (sign < 0.5) {
            sign = -1;

          } else {
            sign = 1;
          }

          if (direction == X)
            x += a * sign;
          else if (direction == Y)
            y += a * sign;
          else
            z += a * sign;
        }

        // Calcola la distanza quadratica dall'origine
        distance = sqrt(x * x + y * y + z * z);
        block_sum += distance;
        
      }

      // Calcola la media nel blocco attuale
      block_sum /= (num_trials / num_blocks);
      block_averages[step - 1] += block_sum;
      block_averages_squared[step - 1] += block_sum * block_sum;
    }

    // Calcola la media tra i blocchi
    block_averages[step - 1] /= num_blocks;
    block_averages_squared[step - 1] /= num_blocks;
    block_errors[step - 1] = sqrt((block_averages_squared[step - 1] - block_averages[step - 1] * block_averages[step - 1]) / num_blocks);
  }

    RandomWalkResult result;
    result.block_averages = block_averages;
    result.block_averages_squared = block_averages_squared;
    result.block_errors = block_errors;
    
    return result;
}

RandomWalkResult random_walk_average_cont(int num_steps, int num_trials, int num_blocks, double a, Random &rnd) {
  vector<double> block_averages(num_steps, 0.0);
  vector<double> block_averages_squared(num_steps, 0.0);
  vector<double> block_errors(num_steps, 0.0);
  double x = 0;
  double y = 0;
  double z = 0;
  double distance_squared = 0;
  double costheta = 0;
  double sintheta = 0;
  double phi = 0;
  double block_sum = 0;

  // Itera su ogni passo del random walk
  for (int step = 1; step <= num_steps; step++) {
    // Itera su ogni blocco
    for (int block = 0; block < num_blocks; block++) {
      block_sum = 0;

      // Itera su ogni lancio nel blocco
      for (int trial = 0; trial < num_trials / num_blocks; trial++) {
        x = 0;
        y = 0;
        z = 0;

        // Esegue il random walk 3D
        for (int s = 0; s < step; s++) {
          phi = rnd.Rannyu(0, 2 * M_PI);
          costheta = rnd.Rannyu(-1,1);
          sintheta = sqrt(1-costheta*costheta);
          x += a * sintheta * cos(phi);
          y += a * sintheta * sin(phi);
          z += a * costheta;
        }

        // Calcola la distanza quadratica dall'origine
        distance_squared = sqrt(x * x + y * y + z * z);
        block_sum += distance_squared;
        
      }

      // Calcola la media nel blocco attuale
      block_sum /= (num_trials / num_blocks);
      block_averages[step - 1] += block_sum;
      block_averages_squared[step - 1] += block_sum * block_sum;
    }

    // Calcola la media tra i blocchi
    block_averages[step - 1] /= num_blocks;
    block_averages_squared[step - 1] /= num_blocks;
    block_errors[step - 1] = sqrt((block_averages_squared[step - 1] - block_averages[step - 1] * block_averages[step - 1]) / num_blocks);
  }

    RandomWalkResult result;
    result.block_averages = block_averages;
    result.block_averages_squared = block_averages_squared;
    result.block_errors = block_errors;
    
    return result;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
