//
//  conftr.c
//
//  Triangular configurations
//
//  Functions:
//    init_neighbor     - Initialize nearest neighbors table
//    write_me          - Magnetization and energy for both spins
//    write_betti       - Betti numbers for positive and negative spins
//    write_euler       - Euler numbers for positive and negative spins
//    write_faces       - Number of k-faces for positive and negative spins
//
//  Copyright (C) 2020 Tak-Shing Chan
//

// Configuration-specific settings
#define NSPIN (L * L)
#define Z 6
#ifndef CELL
#define L2 (L * 2)
#endif

int *neighbor(int i)
{
    static int nn[Z];
    int x = i % L;
    int y = i / L;
    int xm1 = (x - 1 + L) % L;
    int xp1 = (x + 1) % L;
    int ym1 = (y - 1 + L) % L;
    int yp1 = (y + 1) % L;

    // Triangular lattice
    nn[0] = xm1 + y * L;
    nn[1] = xp1 + y * L;
    nn[2] = x + ym1 * L;
    nn[3] = x + yp1 * L;
    nn[4] = xp1 + ym1 * L;
    nn[5] = xm1 + yp1 * L;
    return nn;
}

int pbceuler(int *conf)
{
    // Euler number with periodic boundary conditions
    int euler = 0;
    static int lut[] = {
#ifndef CELL
        0, 1, 1, 0, 1, 0, -2, -1, 1, -2, 0, -1, 0, -1, -1, 0
#else
        0, 1, 2, 0, 2, 0, -2, -1, 1, 2, 0, -2, 0, -2, -1, 0
#endif
    };
    for (int i = 0; i < L; i++) {
        int ip1 = (i + 1) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            euler += lut[conf[i + j * L] * 1 +
                         conf[ip1 + j * L] * 2 +
                         conf[i + jp1 * L] * 4 +
                         conf[ip1 + jp1 * L] * 8];
        }
    }
#ifndef CELL
    return euler / 4;
#else
    return euler / 6;
#endif
}

void write_me(int *conf, FILE *stream)
{
    // Magnetization and energy with periodic boundary conditions
    int M = 0;
    int E = 0;
    for (int i = 0; i < L; i++) {
        int im1 = (i - 1 + L) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            M += conf[i + j * L];
            E -= conf[i + j * L] * (conf[im1 + j * L] +
                                    conf[i + jp1 * L] +
                                    conf[im1 + jp1 * L]);
        }
    }
    fprintf(stream, "%d,%d\n", M, E);
    fflush(stream);
}

void write_betti(int *conf, FILE *stream)
{
#ifndef CELL
    // Betti numbers with periodic boundary conditions
    int b0 = 0;
    int b2 = 1;

    // Depth-first search for connected components
    static int data[NSPIN];
    static int stack[NSPIN];
    memcpy(data, conf, NSPIN * sizeof *conf);
    for (int i = 0; i < NSPIN; i++) {
        switch (data[i]) {
        case 0:
            b2 = 0;             // Not a torus
            break;
        case 1:
            b0++;               // Found connected component
            int sp = 0;         // Initialize stack
            stack[sp++] = i;    // Add spin to stack
            data[i] = 2;        // Mark as visited
            while (sp) {        // Recurse into neighbors
                int spin = stack[--sp];
                int x = spin % L;
                int y = spin / L;
                for (int dx = -1; dx < 2; dx++) {
                    for (int dy = -1; dy < 2; dy++) {
                        int nn = ((x + dx + L) % L) +
                                 ((y + dy + L) % L) * L;
                        if (data[nn] == 1) {
                            stack[sp++] = nn;
                            data[nn] = 2;
                        }
                    }
                }
            }
        }
    }
    fprintf(stream, "%d,%d,%d\n", b0, b0 + b2 - pbceuler(conf), b2);
    fflush(stream);
#endif
}

void write_euler(int *conf, FILE *stream)
{
    // Euler number with periodic boundary conditions
    fprintf(stream, "%d\n", pbceuler(conf));
    fflush(stream);
}

void write_faces(int *conf, FILE *stream)
{
    // Number of k-faces with periodic boundary conditions
    int F0 = 0;
    int F1 = 0;
    int F2 = 0;

#ifndef CELL
    // Create CubeMap with periodic boundary conditions
    static unsigned char cubemap[NSPIN * 4];
    memset(cubemap, 0, NSPIN * 4);
    for (int i = 0; i < L; i++) {
        int i2 = i * 2;
        for (int j = 0; j < L; j++) {
            int j2 = j * 2;
            if (conf[i + j * L]) {
                for (int di = 0; di < 3; di++) {
                    for (int dj = 0; dj < 3; dj++) {
                        cubemap[((i2 + di) % L2) +
                                ((j2 + dj) % L2) * L2] = 1;
                    }
                }
            }
        }
    }

    // Get number of k-faces from CubeMap
    for (int i = 0; i < L; i++) {
        int i2 = i * 2;
        int i2p1 = i2 + 1;
        for (int j = 0; j < L; j++) {
            int j2 = j * 2;
            int j2p1 = j2 + 1;
            F0 += cubemap[i2 + j2 * L2];
            F1 += cubemap[i2 + j2p1 * L2] + cubemap[i2p1 + j2 * L2];
            F2 += cubemap[i2p1 + j2p1 * L2];
        }
    }
#else
    for (int i = 0; i < L; i++) {
        int ip1 = (i + 1) % L;
        for (int j = 0; j < L; j++) {
            int jp1 = (j + 1) % L;
            int quad = conf[i + j * L] * 1 +
                       conf[ip1 + j * L] * 2 +
                       conf[i + jp1 * L] * 4 +
                       conf[ip1 + jp1 * L] * 8;
            F0 += ((quad & 2) == 2);
            F1 += ((quad & 3) == 3) +
                  ((quad & 6) == 6) +
                  ((quad & 10) == 10);
            F2 += ((quad & 7) == 7) +
                  ((quad & 14) == 14);
        }
    }
#endif
    fprintf(stream, "%d,%d,%d\n", F0, F1, F2);
    fflush(stream);
}
