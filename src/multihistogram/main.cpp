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

double avg(vector<double>* d)
{
    vector<double>::iterator cur = d->begin();
    vector<double>::iterator last = d->end();

    double D = 0.0;

    for (; cur != last; cur++)
    {
        D += *cur;
    }

    D /= static_cast<double>(d->size());
    return D;
}

int main(int argc, char** argv)
{
    if (argc != 10)
    {
        cerr << "Usage:\n"
             << "  " << argv[0]
             << " <temperatures_file>"
             << " <L>"
             << " <dimension>"
             << " <results_file>"
             << " <T_min>"
             << " <T_max>"
             << " <n_points>"
             << " <output_file>"
             << " <input_data_dir>\n";

        return 1;
    }

    int seed_rng = rand();
    rlxs_init(1, seed_rng);

    ifstream input(argv[1]);

    if (!input)
    {
        cerr << "Error: cannot open temperature file: "
             << argv[1] << endl;
        return 1;
    }

    const int Ls = atoi(argv[2]);
    const int DIM = atoi(argv[3]);

    ofstream outf(argv[4], ios::out);
    ofstream outf0(argv[8], ios::out);

    if (!outf)
    {
        cerr << "Error: cannot open output file: "
             << argv[4] << endl;
        return 1;
    }

    if (!outf0)
    {
        cerr << "Error: cannot open output file: "
             << argv[8] << endl;
        return 1;
    }

    const double T_min = atof(argv[5]);
    const double T_max = atof(argv[6]);
    const int n = atoi(argv[7]);

    const char* data_dir = argv[9];

    cout << "Tmin = " << T_min
         << ", Tmax = " << T_max << endl;

    int prova;
    double vol = pow(Ls, DIM);

    printf("// Files used:\n");

    string line;
    int nrun = 0;

    while (getline(input, line))
        nrun++;

    input.clear();
    input.seekg(0, input.beg);

    vector<double> Action[nrun];
    vector<double> O[nrun], O1[nrun], O2[nrun];

    double T[nrun];
    string data;

    char fname[512];
    char folder_path[512];

    sprintf(
        folder_path,
        "%s/size_%d",
        data_dir,
        Ls
    );

    for (int i = 0; i < nrun; i++)
    {
        getline(input, data);
        sscanf(data.c_str(), "%lf %d", T + i, &prova);

        sprintf(
            fname,
            "%s/conf-%1.4f-%d.dat",
            folder_path,
            T[i],
            prova
        );

        cout << "Opening file: " << fname << endl;

        ifstream in(fname);

        if (!in)
        {
            cerr << "Error: cannot open file: "
                 << fname << endl;
            return 1;
        }

        double pE, pM;
        double Chip, Chin;
        double F0p, F1p, F2p;
        double F0n, F1n, F2n;

        while (
            in >> pM >> pE >> Chip >> Chin
               >> F0p >> F1p >> F2p
               >> F0n >> F1n >> F2n
        )
        {
            Action[i].push_back(pE);

            O[i].push_back(
                abs(Chip - Chin) / vol
            );

            O1[i].push_back(
                abs(pM) / vol
            );

            O2[i].push_back(
                abs(F1p - F1n) / vol
            );
        }
    }

    double BETA[nrun];

    int nrep = 200;

    // Bootstrap estimates at the original simulation temperatures
    for (int rep = 0; rep < nrep; rep++)
    {
        vector<double> Act_bs[nrun];

        vector<double> O_bs[nrun];
        vector<double> O1_bs[nrun];
        vector<double> O2_bs[nrun];

        vector<double> Osq_bs[nrun];
        vector<double> Osq1_bs[nrun];
        vector<double> Osq2_bs[nrun];

        for (int i = 0; i < nrun; i++)
        {
            int len = Action[i].size();

            double dt[len];
            ranlxs(dt, len);

            double autocorr;
            double autocorr_err;

            BETA[i] = 1.0 / T[i];

            AutoCorr(
                O[i],
                autocorr,
                autocorr_err
            );

            if (autocorr == 0.0)
            {
                cerr << "Warning: autocorrelation time "
                     << "could not be determined at T = "
                     << T[i] << endl;
            }
            else
            {
                cout << setprecision(11)
                     << "Autocorrelation time at T = "
                     << T[i]
                     << " : "
                     << autocorr
                     << endl;
            }

            if (autocorr < 1.0)
                autocorr = 1.0;

            for (int j = 0; autocorr * j < len; j++)
            {
                int idx =
                    static_cast<int>(len * dt[j]);

                Act_bs[i].push_back(
                    Action[i][idx]
                );

                O_bs[i].push_back(
                    O[i][idx]
                );

                Osq_bs[i].push_back(
                    O[i][idx] * O[i][idx]
                );

                O1_bs[i].push_back(
                    O1[i][idx]
                );

                Osq1_bs[i].push_back(
                    O1[i][idx] * O1[i][idx]
                );

                O2_bs[i].push_back(
                    O2[i][idx]
                );

                Osq2_bs[i].push_back(
                    O2[i][idx] * O2[i][idx]
                );
            }
        }

        double Oavr, Osusc;
        double O1avr, O1susc;
        double O2avr, O2susc;

        for (int i = 0; i < nrun; i++)
        {
            Oavr = avg(O_bs + i);

            Osusc =
                avg(Osq_bs + i)
                - Oavr * Oavr;

            O1avr =
                avg(O1_bs + i);

            O1susc =
                avg(Osq1_bs + i)
                - O1avr * O1avr;

            O2avr =
                avg(O2_bs + i);

            O2susc =
                avg(Osq2_bs + i)
                - O2avr * O2avr;

            outf0
                << setprecision(10)
                << rep << " "
                << T[i] << " "
                << Oavr << " "
                << Osusc * vol << " "
                << O1avr << " "
                << O1susc * vol << " "
                << O2avr << " "
                << O2susc * vol << " "
                << endl;
        }
    }

    outf0.close();

    double LogZ[nrun];

    for (int i = 0; i < nrun; i++)
    {
        LogZ[i] = 0.0;
    }

    // Multihistogram reweighting
    for (int rep = 0; rep < nrep; rep++)
    {
        vector<double> Act_bs[nrun];

        vector<double> O_bs[nrun];
        vector<double> O1_bs[nrun];
        vector<double> O2_bs[nrun];

        for (int i = 0; i < nrun; i++)
        {
            int len = Action[i].size();

            double dt[len];
            ranlxs(dt, len);

            double autocorr;
            double autocorr_err;

            AutoCorr(
                O[i],
                autocorr,
                autocorr_err
            );

            if (autocorr == 0.0)
            {
                cerr << "Warning: autocorrelation time "
                     << "could not be determined at beta = "
                     << BETA[i] << endl;
            }
            else
            {
                cout << setprecision(11)
                     << "Autocorrelation time at beta = "
                     << BETA[i]
                     << " : "
                     << autocorr
                     << endl;
            }

            if (autocorr < 1.0)
                autocorr = 1.0;

            for (int j = 0; autocorr * j < len; j++)
            {
                int idx =
                    static_cast<int>(len * dt[j]);

                Act_bs[i].push_back(
                    Action[i][idx]
                );

                O_bs[i].push_back(
                    O[i][idx]
                );

                O1_bs[i].push_back(
                    O1[i][idx]
                );

                O2_bs[i].push_back(
                    O2[i][idx]
                );
            }
        }

        MultiHistRw(
            nrun,
            BETA,
            Act_bs,
            LogZ
        );

        double b;
        double FreeNRG;

        double Oavr, Osusc;
        double O1avr, O1susc;
        double O2avr, O2susc;

        // beta decreases when temperature increases
        const double B_min = 1.0 / T_max;
        const double B_max = 1.0 / T_min;

        for (int i = 0; i < n; i++)
        {
            b =
                ((B_max - B_min)
                 / static_cast<double>(n - 1))
                * i
                + B_min;

            FreeNRG =
                LnZ(
                    b,
                    nrun,
                    BETA,
                    Act_bs,
                    LogZ
                );

            Oavr =
                exp(
                    LnO_n(
                        O_bs,
                        1.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                );

            Osusc =
                exp(
                    LnO_n(
                        O_bs,
                        2.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                )
                - Oavr * Oavr;

            O1avr =
                exp(
                    LnO_n(
                        O1_bs,
                        1.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                );

            O1susc =
                exp(
                    LnO_n(
                        O1_bs,
                        2.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                )
                - O1avr * O1avr;

            O2avr =
                exp(
                    LnO_n(
                        O2_bs,
                        1.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                );

            O2susc =
                exp(
                    LnO_n(
                        O2_bs,
                        2.0,
                        b,
                        nrun,
                        BETA,
                        Act_bs,
                        LogZ
                    )
                    - FreeNRG
                )
                - O2avr * O2avr;

            outf
                << setprecision(10)
                << rep << " "
                << 1.0 / b << " "
                << FreeNRG << " "
                << Oavr << " "
                << Osusc * vol << " "
                << O1avr << " "
                << O1susc * vol << " "
                << O2avr << " "
                << O2susc * vol << " "
                << endl;
        }
    }

    outf.close();

    return 0;
}