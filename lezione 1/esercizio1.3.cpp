#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

int main() {
    Random rnd;

    double d = 10;    // distanza tra le righe
    double lungh = 5; // lunghezza ago
    double throws = 1000000;
    double nblocchi = 100;
    double counter = 0;
    double xi, xf = 0;
    double dx;
    double L = throws/nblocchi;


    /*for (int i = 0; i < throws; i++) {
        xc = rnd.Rannyu(0, d);
        theta = rnd.Rannyu(-M_PI / 2, M_PI / 2);
        xf = xc + (lungh / 2) * cos(theta);
        xi = xc - (lungh / 2) * cos(theta);
        
        if (xi < 0 || xf > d) {
            counter++;
        }
    }*/

    vector<double> pi (nblocchi);
    
   for(int i=0; i<nblocchi; i++) {

        xi =0 ;
        dx = 0;
        xf = 0;

        for(int j=0; j<L; j++) {
        xi = rnd.Rannyu(0, d);
        dx = lungh*rnd.Cos();

        xf = xi+dx;

            if(xf<0 || xf>d) {
    
                counter ++;
            }

        }
        
        pi[i] = ((2 * lungh * L) / (counter * d));
        counter = 0;

    }

    vector<double> pi_mediablocchi = mediablocchi(nblocchi,throws,pi);
    vector<double> pi_mediablocchi2 = mediablocchi2(nblocchi,throws,pi);
    vector<double> pi_err = errore(nblocchi,pi_mediablocchi, pi_mediablocchi2);

   fstream foutput; //foutput è il nome che do alla variabile
   foutput.open("pi.dat", ios::out); //risultato.dat è il nome del file
   if (foutput.good()) //controllo che funzioni e gli faccio stampare, altrimenti gli faccio dare errore
    {
      for(double i=0; i<nblocchi; i++) {
      foutput << pi_mediablocchi[i] << " " << pi_err[i] << endl;
      }
    }
  else
    cout << "Errore: file is not good!" << endl;

  foutput.close(); //chiudo il file

   /*for(int i=0; i<throws; i++) {
        xi = rnd.Rannyu(0, d);
        dx = lungh*rnd.Cos();

        xf = xi+dx;

            if(xf<0 || xf>d) {
    
                counter ++;
            }

    }

    double stima_pi = ((2 * lungh * throws) / (counter * d));

    cout << "Stima di pi: " << stima_pi << endl;*/

    return 0;
}



