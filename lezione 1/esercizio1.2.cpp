#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

using namespace std;

int main() {

    Random rnd; //Creazione oggetto random
    int throws = 10000; //Numero di lanci

    // Generiamo i dati per diversi valori di N e li scriviamo nei file corrispondenti
    generateData("n1.dat", throws, 1, rnd);
    generateData("n2.dat", throws, 2, rnd);
    generateData("n10.dat", throws, 10, rnd);
    generateData("n100.dat", throws, 100, rnd);

    return 0;
}
