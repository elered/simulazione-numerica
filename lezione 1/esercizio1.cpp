#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;
 
int main (){

    Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    int M = 100000;              // Total number of throws
    int nblocchi = 100;                 // Number of blocks
    double L = M/nblocchi;           //Number of throws in each block, please use for M a multiple of N
    double mediablocco = 0;
    double mediablocco2 = 0;
    double num = 0;
    vector <double> var(nblocchi);
    vector <double> Var(nblocchi);
    double sum = 0;         //somme parziali
    double sum2 = 0;
    vector <double> vec_media(nblocchi);
    vector <double> vec_media2(nblocchi);
    vector <double> sumprog(nblocchi);
    vector <double> sumprog2(nblocchi);
    vector <double> Sumprog(nblocchi);
    vector <double> Sumprog2(nblocchi);
    vector <double> ave(nblocchi);
    vector <double> ave2(nblocchi);

   rnd.SaveSeed();
    
    for(int k=0; k<nblocchi; k++) {
      sum = 0;
      sum2 = 0;
      mediablocco2=0;
      mediablocco=0;

        for(int i=0; i<L; i++) {
        num = rnd.Rannyu();
        sum += num;
        sum2 += pow(num -0.5,2);
        }
        mediablocco = sum/L;
        mediablocco2 = mediablocco*mediablocco;
        vec_media[k]=mediablocco;
        vec_media2[k]= mediablocco2;      
        ave[k]=sum2/L;
        ave2[k]=ave[k]*ave[k];
    }
  
   for(int i=0; i<nblocchi; i++) {
      for (int j=0; j<=i; j++) {

         sumprog[i] += vec_media[j];
         sumprog2[i] += vec_media2[j];
      }
      sumprog[i]/=(i+1);
      sumprog2[i]/=(i+1);
      if(i==0){
         var[i]=0;
      } else {
      var[i] = sqrt((sumprog2[i] - sumprog[i] * sumprog[i]) / (i));
      }
   }

   for(int i=0; i<nblocchi; i++) {
      for (int j=0; j<=i; j++) {

         Sumprog[i] += ave[j];
         Sumprog2[i] += ave2[j];
      }
      Sumprog[i]/=(i+1);
      Sumprog2[i]/=(i+1);
      if(i==0){
         Var[i]=0;
      } else {
      Var[i] = sqrt((Sumprog2[i] - Sumprog[i] * Sumprog[i]) / (i));
      }
   }

   fstream foutput; //foutput è il nome che do alla variabile
   foutput.open("risultato.dat", ios::out); //risultato.dat è il nome del file
   if (foutput.good()) //controllo che funzioni e gli faccio stampare, altrimenti gli faccio dare errore
    {
      for(int i=0; i<nblocchi; i++) {
      foutput << sumprog[i] << " " << var[i] << " " << Sumprog[i] << " " << Var[i] << endl;
      }
    }
  else
    cout << "Errore: file is not good!" << endl;

  foutput.close(); //chiudo il file


    return 0;
}

   