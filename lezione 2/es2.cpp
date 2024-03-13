#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

int main() {
  // Inizializzazione oggetto Random per la generazione di numeri casuali
  Random rnd;

  // Definizione dei parametri della simulazione
  double throws = 10000;
  double passi = 100;
  double blocchi = 100;
  double a = 1;

  RandomWalkResult result = random_walk_average(passi, throws, blocchi, a, rnd);
    vector<double> block_averages_discr = result.block_averages;
    vector<double> block_averages_discr_squared = result.block_averages_squared;
    vector<double> block_errors_discr = result.block_errors;

    RandomWalkResult result_cont = random_walk_average_cont(passi, throws, blocchi, a, rnd);
    vector<double> block_averages_cont = result.block_averages;
    vector<double> block_averages_cont_squared = result.block_averages_squared;
    vector<double> block_errors_cont = result.block_errors;


  // Stampa dei dati su file
  stampa("dati.dat", blocchi, block_averages_discr, block_errors_discr, block_averages_cont, block_errors_cont);

  return 0;
}