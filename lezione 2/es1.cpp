#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "random.h"

using namespace std;

int main() {
  Random rnd;
  double max = 1;
  double min = 0;
  double throws = 10000;
  double blocchi = 100;
  vector<double> ave_blocchi = integraleave(min, max, throws, blocchi, rnd);
  vector<double> mediablocchi_int = mediablocchi(blocchi, throws, ave_blocchi);
  vector<double> mediablocchi2_int =
      mediablocchi2(blocchi, throws, ave_blocchi);
  vector<double> errore_int =
      errore(blocchi, mediablocchi_int, mediablocchi2_int);

  vector<double> avehom_blocchi =
      integralehom(min, max, 3 / 2, rnd, throws, blocchi);
  vector<double> avehomblocchi_int =
      mediablocchi(blocchi, throws, avehom_blocchi);
  vector<double> avehomblocchi2_int =
      mediablocchi2(blocchi, throws, avehom_blocchi);
  vector<double> errorehom_int =
      errore(blocchi, avehomblocchi_int, avehomblocchi2_int);

  stampa("risultati.dat", blocchi, mediablocchi_int, errore_int,
         avehomblocchi_int, errorehom_int);

  return 0;
}