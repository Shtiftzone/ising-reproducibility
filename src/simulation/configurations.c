#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ranlux.h"

#define NRAND 360

#define L 512               // Rozmiar kratki LxL
#define NSPIN (L * L)      // Liczba spinów
#define Z 4                // Liczba s¹siadów (2D kratka)
#define NTHERM 10        // Liczba kroków termalizacji

int s[NSPIN];             // Konfiguracja spinów
int stack[NSPIN];         // Stos dla algorytmu Wolffa

// Funkcja generuj¹ca s¹siadów
int *neighbor(int i) {
    static int nn[Z];
    int x = i % L;
    int y = i / L;
    nn[0] = x + ((y + 1) % L) * L; // S¹siad góra
    nn[1] = x + ((y - 1 + L) % L) * L; // S¹siad dó³
    nn[2] = ((x + 1) % L) + y * L; // S¹siad prawo
    nn[3] = ((x - 1 + L) % L) + y * L; // S¹siad lewo
    return nn;
}

void sranlux(const char *filename)
{
    // Get seed from file
    unsigned seed;
    FILE *stream = fopen(filename, "r");
    if (stream && fscanf(stream, "%u\n", &seed) == 1) {
        fclose(stream);
        rlxd_init(1, seed);
    } else {
        // If failed, create a new file
        stream = fopen("/dev/urandom", "r");
        if (stream == NULL || fread(&seed, sizeof seed, 1, stream) != 1) {
            fprintf(stderr, "Couldn't read from /dev/urandom\n");
            exit(EXIT_FAILURE);
        }
        fclose(stream);
        seed &= 0x3fffffff;
        rlxd_init(1, seed);
        stream = fopen(filename, "w");
        if (stream == NULL || fprintf(stream, "%u\n", seed) < 0) {
            fprintf(stderr, "Couldn't write to %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fclose(stream);
    }
}

double dranlux(void)
{
    // LIFO buffering of random numbers with skip
    static int irand = 0;
    static double ran[NRAND];
    if (irand == 0) {
        ranlxd(ran, NRAND);
        irand = NRAND - 1;
    }
    return ran[irand--];
}


void metro_update(double T) {
    for (int i = 0; i < NSPIN; i++) {
        int *nn = neighbor(i);
        int sum = 0;
        for (int j = 0; j < Z; j++)
            sum += s[nn[j]];
        int dE = 2 * s[i] * sum;
        if (dE < 0 || dranlux() < exp(-dE / T))
            s[i] = -s[i];
    }
}

void wolff_update(double p) {
    int seed = (int)(dranlux() * NSPIN);
    int sp = 0;
    stack[sp++] = seed;
    s[seed] = -s[seed];
    while (sp) {
        int site = stack[--sp];
        int *nn = neighbor(site);
        for (int j = 0; j < Z; j++) {
            if (s[site] == -s[nn[j]] && dranlux() < p) {
                stack[sp++] = nn[j];
                s[nn[j]] = s[site];
            }
        }
    }
}

void save_configuration(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Couldn't open file %s for writing\n", filename);
        exit(EXIT_FAILURE);
    }
    for (int y = 0; y < L; y++) {
        for (int x = 0; x < L; x++) {
            fprintf(file, "%d ", s[x + y * L]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

int main() {
    double T_start = 2.0;
    double T_end = 2.5;
    double dT = 0.005;
    int config_number = 0;

    rlxd_init(1, 12345);  // Inicjalizacja generatora losowego

    // Inicjalizacja konfiguracji spinów
    for (int i = 0; i < NSPIN; i++)
        s[i] = (dranlux() < 0.5) ? 1 : -1;

    for (double T = T_start; T <= T_end; T += dT) {
        double p = 1.0 - exp(-2.0 / T);

        // Termalizacja
        for (int i = 0; i < NTHERM; i++) {
            metro_update(T);
            wolff_update(p);

            // Zapis konfiguracji po ka¿dej iteracji
            char filename[100];
            sprintf(filename, "/users/project1/pt01192/configurations/config_T%.3f_iter_%03d.txt", T, i);
            save_configuration(filename);
        }

    }

    return 0;
}