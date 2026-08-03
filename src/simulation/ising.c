//
//  ising.c
//
//  Ising model with hybrid Metropolis and Wolff updates
//
//  Parameters:
//    T                 - Temperature
//    Nconf             - Number of configurations
//
//  Inputs:
//    Sfile             - Seed
//    Ndecorr           - Number of decorrelation steps (optional)
//
//  Outputs:
//    Mefile            - Magnetization and energy for both spins
//    Bpfile, Bnfile    - Betti numbers for positive and negative spins
//    Epfile, Enfile    - Euler numbers for positive and negative spins
//    Fpfile, Fnfile    - Number of k-faces for positive and negative spins
//
//  Copyright (C) 2020 Tak-Shing Chan and Davide Vadacchino
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ranlux.h"

#define NRAND 360

#ifndef DIM
#include "conftr.c"
#elif DIM == 2
#include "conf2d.c"
#elif DIM == 3
#include "conf3d.c"
#elif DIM == 4
#include "conf4d.c"
#else
#error Unsupported configuration.
#endif

int *s;         // Spin configuration
int *stack;     // Pseudo-stack area

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

FILE *file_notna(const char *filename)
{
    // Drop missing values
    return strcmp(filename, "NA") ? fopen(filename, "w") : NULL;
}

#ifdef DIM
int *neighbor(int i)
{
    // Von Neumann neighborhood on a cubic lattice
    static int nn[Z];
    int p = 1;
    int q = 1 - L;
    int r = i;
    for (int j = 0; j < Z; j += 2) {
        nn[j] = (r + 1) % L ? i + p : i + q;
        nn[j + 1] = r % L ? i - p : i - q;
        p *= L;
        q *= L;
        r /= L;
    }
    return nn;
}
#endif

void metro_update(double T)
{
    // Sequential Metropolis updates
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

void wolff_update(int N, double p)
{
    // Random Wolff updates
    for (int i = 0; i < N; i++) {
        int seed = (int) (dranlux() * NSPIN);
        int sp = 0;             // Initialize stack
        stack[sp++] = seed;     // Add seed to stack
        s[seed] = -s[seed];     // Reverse spin
        while (sp) {
            seed = stack[--sp]; // Remove spin from stack
            int *nn = neighbor(seed);   // Cycle over neighbors
            for (int j = 0; j < Z; j++) {
                if (dranlux() < p && s[seed] == -s[nn[j]]) {
                    stack[sp++] = nn[j];
                    s[nn[j]] = s[seed];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    // Empirically-determined constants
    int Ntherm = 2000;
    int Nwolff = 15;
    int Ndecorr = 4;
    if (argc < 11) {
        fprintf(stderr,
                "%s T Nconf Sfile Mefile Bpfile Bnfile Epfile Enfile Fpfile Fnfile [Ndecorr]\n",
                argv[0]);
        return EXIT_FAILURE;
    }
    s = malloc(NSPIN * sizeof *s);
    stack = malloc(NSPIN * sizeof *stack);
    if (!s || !stack) {
        fprintf(stderr, "Couldn't allocate memory\n");
        return EXIT_FAILURE;
    }
    double T = strtod(argv[1], NULL);
    double p = 1.0 - exp(-2.0 / T);
    int Nconf = (int) strtol(argv[2], NULL, 10);
    sranlux(argv[3]);
    FILE *Mefile = file_notna(argv[4]);
    FILE *Bpfile = file_notna(argv[5]);
    FILE *Bnfile = file_notna(argv[6]);
    FILE *Epfile = file_notna(argv[7]);
    FILE *Enfile = file_notna(argv[8]);
    FILE *Fpfile = file_notna(argv[9]);
    FILE *Fnfile = file_notna(argv[10]);
    if (argc == 12)
        Ndecorr = (int) strtol(argv[11], NULL, 10);

    // Initialization and thermalization
    for (int i = 0; i < NSPIN; i++)
        s[i] = dranlux() < 0.5 ? 1 : -1;
    for (int i = 0; i < Ntherm; i++) {
		metro_update(T);
		wolff_update(Nwolff, p);

		//Save configuration to file
		/*#ifdef DIM
		if (DIM == 2) {
			char filename[50];
			sprintf(filename, "spin_config_%04d.txt", i);
			FILE *spinfile = fopen(filename, "w");
			for (int x = 0; x < L; x++) {
				for (int y = 0; y < L; y++) {
					fprintf(spinfile, "%d ", s[x + y * L]);
				}
				fprintf(spinfile, "\n");
			}
			fclose(spinfile);
		}
		#endif*/
}

    // Starting simulation
    while (Nconf--) {
        // For both spins
        if (Mefile)
            write_me(s, Mefile);

        // For positive spins
        for (int i = 0; i < NSPIN; i++)
            stack[i] = s[i] > 0;
        if (Bpfile)
            write_betti(stack, Bpfile);
        if (Epfile)
            write_euler(stack, Epfile);
        if (Fpfile)
            write_faces(stack, Fpfile);

        // For negative spins
        for (int i = 0; i < NSPIN; i++)
            stack[i] = s[i] < 0;
        if (Bnfile)
            write_betti(stack, Bnfile);
        if (Enfile)
            write_euler(stack, Enfile);
        if (Fnfile)
            write_faces(stack, Fnfile);

        // Decorrelation
        for (int i = 0; i < Ndecorr; i++) {
            metro_update(T);
            wolff_update(Nwolff, p);
        }
    }
    return 0;
}
