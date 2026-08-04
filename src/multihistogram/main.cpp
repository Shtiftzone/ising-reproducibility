#include "ranlxs.h"
#include "MultiHistRW.h"
#include "Autocorr.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

double avg(vector<double> * d)
{
    vector<double>::iterator cur = d->begin();
    vector<double>::iterator last = d->end();
    double D = 0.;
    for(; cur != last; cur++)
    {
        D += *cur;
    }
    D /= (double)d->size();
    return D;
}

int main(int argc, char** argv) {
    int seed_rng = rand();
    rlxs_init(1, seed_rng);

    ifstream input(argv[1]);
    const int Ls = atoi(argv[2]); // Rozmiar siatki (size_L)
    const int DIM = atoi(argv[3]);

    ofstream outf;
    outf.open(argv[4], ios::out);
    ofstream outf0;
    outf0.open(argv[8], ios::out);

    const double T_min = atof(argv[6]);
    const double T_max = atof(argv[5]);
    cout << "Tmin = " << T_min << ", Tmax = " << T_max << endl;
    const int n = atoi(argv[7]); // Iloœæ punktów

    int prova;
    double vol = pow(Ls, DIM);
    printf("//Files used:\n");

    string line;
    int nrun = 0;
    while (getline(input, line)) nrun++;
    input.clear();
    input.seekg(0, input.beg);

    vector<double> Action[nrun];
    vector<double> O[nrun], O1[nrun], O2[nrun];

    double T[nrun];
    string data;
    char fname[512]; // Wiêkszy rozmiar dla pe³nej œcie¿ki
    double tmpB;

    // Tworzymy poprawn¹ œcie¿kê do folderu `size_L`
    char folder_path[256];
    sprintf(folder_path, "data/eul2d/size_%d", Ls); // Tworzymy œcie¿kê dla odpowiedniego rozmiaru

    for (int i = 0; i < nrun; i++) {
        getline(input, data);
        sscanf(data.c_str(), "%lf %d\n", T + i, &prova);

        // Konstruowanie poprawnej œcie¿ki pliku dla `size_L`
        sprintf(fname, "%s/conf-%1.4f-%d.dat", folder_path, T[i], prova);

        cout << "Otwieram plik: " << fname << endl;
        ifstream in(fname);
        if (!in) {
            cerr << "B³¹d: Nie mo¿na otworzyæ pliku " << fname << endl;
            exit(1);
        }

        double pE, pM;
        double Chip, Chin;
        double F0p, F1p, F2p;
        double F0n, F1n, F2n;

        while (in >> pM >> pE >> Chip >> Chin >>
                      F0p >> F1p >> F2p >>
                      F0n >> F1n >> F2n) {
            Action[i].push_back(pE);
            O[i].push_back(abs(Chip - Chin) / vol);
            O1[i].push_back(abs(pM) / vol);
            O2[i].push_back(abs(F1p - F1n) / vol);
        }
    }

  double BETA[nrun];
  int nrep=200; //// of bootstrap samples
  for(int rep=0; rep<nrep; rep++) {
    vector<double> Act_bs[nrun];
    vector<double> O_bs[nrun];
    vector<double> O1_bs[nrun];
    vector<double> O2_bs[nrun];
    vector<double> Osq_bs[nrun];
    vector<double> Osq1_bs[nrun];
    vector<double> Osq2_bs[nrun];
  
    //Create Bootstrap sample
    for (int i=0; i<nrun; i++) {
      int len = Action[i].size();
      double dt[len];
      ranlxs(dt, len);

      double autocorr, autocorr_err;
      BETA[i] = 1./T[i];
      AutoCorr(O[i], autocorr, autocorr_err);
      //cout << " autocorr time = " << autocorr << endl;
      if(autocorr==0.) //Non si riesce a determinare il tempo di autocorrelazione
        cerr<<"ATTENZIONE: autocorr -> "<< T[i]<<endl;
      else cout << setprecision(11) << "Autocorr at beta = " << T[i] << " : " << autocorr << endl;
      if(autocorr<1.) autocorr=1.0;
    
      for(int j=0; autocorr*j<len; j++) {
        int idx = static_cast<int>(len*dt[j]);
        Act_bs[i].push_back(Action[i][idx]);
        O_bs[i].push_back(O[i][idx]);
        Osq_bs[i].push_back(O[i][idx]*O[i][idx]);

        O1_bs[i].push_back(O1[i][idx]);
        Osq1_bs[i].push_back(O1[i][idx]*O1[i][idx]);

        O2_bs[i].push_back(O2[i][idx]);
        Osq2_bs[i].push_back(O2[i][idx]*O2[i][idx]);
        //cout <<fixed<< " E_bs = " << Act_bs[i][j] << " , O_bs = " << O_bs[i][j] << endl; 
      }
    }

    double Oavr, Osusc;
    double O1avr, O1susc;
    double O2avr, O2susc;
    for(int i=0; i<nrun; i++) {
      Oavr  =  avg( O_bs+i );
      Osusc =  avg( Osq_bs+i ) - Oavr*Oavr;

      O1avr =  avg( O1_bs+i );
      O1susc=  avg( Osq1_bs+i ) - O1avr*O1avr;
      
      O2avr =  avg( O2_bs+i );
      O2susc=  avg( Osq2_bs+i ) - O2avr*O2avr;

      outf0<<setprecision(10)<<rep<<" "
        << T[i] <<" "
        <<Oavr<<" "
        <<Osusc*vol<<" "
        <<O1avr<<" "
        <<O1susc*vol<<" "
        <<O2avr<<" "
	<<O2susc*vol<<" "
        <<endl;
    }
  }
  outf0.close();


  double LogZ[nrun];
  for(int i=0; i<nrun; i++) {
    LogZ[i]=0.0; //Initial guess for Z
  }

  for(int rep=0; rep<nrep; rep++) {
    vector<double> Act_bs[nrun];
    vector<double> O_bs[nrun];
    vector<double> O1_bs[nrun];
    vector<double> O2_bs[nrun];
    
    //Create Bootstrap sample
    for (int i=0; i<nrun; i++) {
      int len = Action[i].size();
      double dt[len];
      ranlxs(dt, len);

      double autocorr, autocorr_err;
      AutoCorr(O[i], autocorr, autocorr_err);
      //cout << " autocorr time = " << autocorr << endl;
      if(autocorr==0.) //Non si riesce a determinare il tempo di autocorrelazione
	cerr<<"ATTENZIONE: autocorr O -> "<< BETA[i]<<endl;
      else cout << setprecision(11) << "Autocorr at beta = " << BETA[i] << " : " << autocorr << endl;
      if(autocorr<1.) autocorr=1.0;
      
      for(int j=0; autocorr*j<len; j++) {
	int idx = static_cast<int>(len*dt[j]);
	Act_bs[i].push_back(Action[i][idx]);
	O_bs[i].push_back(O[i][idx]);
	O1_bs[i].push_back(O1[i][idx]);
	O2_bs[i].push_back(O2[i][idx]);
	//cout <<fixed<< " E_bs = " << Act_bs[i][j] << " , O_bs = " << O_bs[i][j] << endl; 
      }
    }

    //Perform Reweigthing
    MultiHistRw(nrun, BETA, Act_bs, LogZ);
    //Calculate average and susceptibility at intermediate betas
    double b, FreeNRG, Oavr, Osusc;
    double O1avr, O1susc;
    double O2avr, O2susc;
    double B_min= 1./T_max;
    double B_max= 1./T_min;
    for(int i=0; i<n; i++) {
      b=((B_max-B_min)/(double)(n-1))*i+B_min;      
      FreeNRG=LnZ(b, nrun, BETA, Act_bs, LogZ);
      Oavr=exp(LnO_n(O_bs, 1.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG);
      Osusc=(exp(LnO_n(O_bs, 2.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG)-Oavr*Oavr);
      O1avr=exp(LnO_n(O1_bs, 1.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG);
      O1susc=(exp(LnO_n(O1_bs, 2.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG)-O1avr*O1avr);
      O2avr=exp(LnO_n(O2_bs, 1.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG);
      O2susc=(exp(LnO_n(O2_bs, 2.0, b, nrun, BETA, Act_bs, LogZ)-FreeNRG)-O2avr*O2avr);

      outf<<setprecision(10)<<rep<<" "
	  << 1./b <<" "
	  <<FreeNRG <<" "

	  <<Oavr <<" "
	  <<Osusc*vol <<" "

	  <<O1avr <<" "
	  <<O1susc*vol <<" "

	  <<O2avr <<" "
	  <<O2susc*vol <<" "
	  <<endl;
    }

  }
  outf.close();

  return 0;

}
