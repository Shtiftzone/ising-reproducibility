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
//  Modifications Copyright (C) 2026 Mateusz Maslowski
//

// Configuration-specific settings
#define NSPIN (L * L)
#define Z 6


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

#ifndef SPIN_AS_VERTEX

        /*
         * Hexagonal cell complex.
         *
         * Each spin is represented as a hexagonal 2-cell.
         *
         * Local block:
         *
         *              2
         *           /     \
         *          1       4
         *           \     /
         *              8
         *
         * Bits:
         *
         *   hex1 -> 1
         *   hex2 -> 2
         *   hex4 -> 4
         *   hex8 -> 8
         *
         * The block is translated over the complete periodic lattice.
         *
         * Multiplicities under this covering:
         *
         *   vertex             : 1
         *   horizontal edge    : 1
         *   diagonal edge      : 2
         *   hexagonal face     : 4
         *
         * Therefore the local integer contribution is
         *
         *   4*V - 4*E_horizontal - 2*E_diagonal + F
         *
         * and the final sum is divided by 4.
         */

         0,  1,  1,  0,
         1,  2,  0, -1,
         1,  0, -2, -1,
         0, -1, -1,  0

#else

        /*
         * Spin-as-vertex simplicial complex
         * on the triangular lattice.
         */

         0,  1,  2,  0,
         2,  0, -2, -1,
         1,  2,  0, -2,
         0, -2, -1,  0

#endif
    };


#ifndef SPIN_AS_VERTEX

    /*
     * Spin-as-cell representation.
     *
     * Every lattice site is used successively as hexagon 2:
     *
     *              2
     *           /     \
     *          1       4
     *           \     /
     *              8
     *
     * Coordinates:
     *
     *   hex1 = (i-1, j+1)
     *   hex2 = (i,   j)
     *   hex4 = (i+1, j)
     *   hex8 = (i,   j+1)
     */

    for (int i = 0; i < L; i++) {

        int im1 = (i - 1 + L) % L;
        int ip1 = (i + 1) % L;

        for (int j = 0; j < L; j++) {

            int jp1 = (j + 1) % L;

            int quad =
                conf[im1 + jp1 * L] * 1 +
                conf[i   + j   * L] * 2 +
                conf[ip1 + j   * L] * 4 +
                conf[i   + jp1 * L] * 8;

            euler += lut[quad];
        }
    }

    return euler / 4;


#else

    /*
     * Spin-as-vertex simplicial construction.
     */

    for (int i = 0; i < L; i++) {

        int ip1 = (i + 1) % L;

        for (int j = 0; j < L; j++) {

            int jp1 = (j + 1) % L;

            euler += lut[
                conf[i   + j   * L] * 1 +
                conf[ip1 + j   * L] * 2 +
                conf[i   + jp1 * L] * 4 +
                conf[ip1 + jp1 * L] * 8
            ];
        }
    }

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

            E -= conf[i + j * L] *
                 (conf[im1 + j   * L] +
                  conf[i   + jp1 * L] +
                  conf[im1 + jp1 * L]);
        }
    }

    fprintf(stream, "%d,%d\n", M, E);
    fflush(stream);
}


void write_betti(int *conf, FILE *stream)
{
#ifndef SPIN_AS_VERTEX

    /*
     * Betti numbers of the spin-as-cell hexagonal complex
     * with periodic boundary conditions.
     *
     * b0 = number of connected components of active hexagons
     *
     * b2 = 1 only if all hexagonal cells are active,
     *      in which case the complete torus is present
     *
     * b1 follows from Euler-Poincare:
     *
     *      chi = b0 - b1 + b2
     *
     * hence
     *
     *      b1 = b0 + b2 - chi
     */

    int b0 = 0;
    int b2 = 1;

    static int data[NSPIN];
    static int stack[NSPIN];

    memcpy(data, conf, NSPIN * sizeof *conf);

    for (int i = 0; i < NSPIN; i++) {

        switch (data[i]) {

        case 0:
            // At least one hexagonal cell is absent,
            // so the complete torus is not filled.
            b2 = 0;
            break;

        case 1:
        {
            // Found a new connected component.
            b0++;

            int sp = 0;

            stack[sp++] = i;
            data[i] = 2;

            while (sp) {

                int spin = stack[--sp];

                /*
                 * Two hexagonal cells share an edge exactly when
                 * their centers are nearest neighbours on the
                 * triangular lattice.
                 */
                int *nn = neighbor(spin);

                for (int k = 0; k < Z; k++) {

                    int next = nn[k];

                    if (data[next] == 1) {

                        data[next] = 2;
                        stack[sp++] = next;
                    }
                }
            }

            break;
        }

        default:
            // Already visited.
            break;
        }
    }

    int chi = pbceuler(conf);
    int b1 = b0 + b2 - chi;

    fprintf(stream, "%d,%d,%d\n", b0, b1, b2);
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

#ifndef SPIN_AS_VERTEX

    /*
     * Spin-as-cell hexagonal complex.
     *
     * Each lattice site corresponds to one hexagonal 2-cell.
     *
     * Every spin is used once as the anchor hexagon 2:
     *
     *              2
     *           /     \
     *          1       4
     *           \     /
     *              8
     *
     * The anchor owns:
     *
     *   one 2-face:
     *       2
     *
     *   three edges:
     *       (1,2)
     *       (2,4)
     *       (2,8)
     *
     *   two vertices:
     *       (1,2,8)
     *       (2,4,8)
     *
     * Translating the anchor over all sites counts every
     * hexagonal face, edge and vertex exactly once.
     *
     * An edge is present if at least one of its two incident
     * hexagons is active.
     *
     * A vertex is present if at least one of its three incident
     * hexagons is active.
     */

    for (int i = 0; i < L; i++) {

        int im1 = (i - 1 + L) % L;
        int ip1 = (i + 1) % L;

        for (int j = 0; j < L; j++) {

            int jp1 = (j + 1) % L;

            int hex1 = conf[im1 + jp1 * L];
            int hex2 = conf[i   + j   * L];
            int hex4 = conf[ip1 + j   * L];
            int hex8 = conf[i   + jp1 * L];


            // One hexagonal 2-face owned by the anchor.
            F2 += hex2;


            // Three edges owned by the anchor.
            F1 += (hex1 || hex2) +
                  (hex2 || hex4) +
                  (hex2 || hex8);


            // Two vertices owned by the anchor.
            F0 += (hex1 || hex2 || hex8) +
                  (hex2 || hex4 || hex8);
        }
    }

#else

    /*
     * Spin-as-vertex simplicial construction.
     *
     * Each active spin is a vertex.
     * Edges and triangular faces are included when all of
     * their vertices are active.
     */

    for (int i = 0; i < L; i++) {

        int ip1 = (i + 1) % L;

        for (int j = 0; j < L; j++) {

            int jp1 = (j + 1) % L;

            int quad =
                conf[i   + j   * L] * 1 +
                conf[ip1 + j   * L] * 2 +
                conf[i   + jp1 * L] * 4 +
                conf[ip1 + jp1 * L] * 8;

            // One anchor vertex.
            F0 += ((quad & 2) == 2);

            // Three edges incident to the anchor.
            F1 += ((quad & 3)  == 3) +
                  ((quad & 6)  == 6) +
                  ((quad & 10) == 10);

            // Two triangular faces.
            F2 += ((quad & 7)  == 7) +
                  ((quad & 14) == 14);
        }
    }

#endif

    fprintf(stream, "%d,%d,%d\n", F0, F1, F2);
    fflush(stream);
}