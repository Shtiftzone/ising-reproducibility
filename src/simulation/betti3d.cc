//
//  betti3d.cc
//
//  3D Betti numbers with periodic boundary conditions
//
//  Copyright (C) 2020 Tak-Shing Chan
//

#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Persistent_cohomology.h>
#include <vector>
#include <algorithm>

using namespace std;

extern "C" {
    void betti3d(int L, int *inp, int *outp)
    {
        typedef Gudhi::cubical_complex::Bitmap_cubical_complex_periodic_boundary_conditions_base<double> Bitmap_base;
        typedef Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_base> Bitmap_cubical_complex;
        typedef Gudhi::persistent_cohomology::Field_Zp Field_Zp;
        typedef Gudhi::persistent_cohomology::Persistent_cohomology<Bitmap_cubical_complex, Field_Zp> Persistent_cohomology;

        vector<unsigned> sizes(3, L);
        vector<double> data(inp, inp + L * L * L);
        vector<bool> periodic(3, true);
        Bitmap_cubical_complex b(sizes, data, periodic);
        Persistent_cohomology pcoh(b, true);
        pcoh.init_coefficients(11);
        pcoh.compute_persistent_cohomology();
        vector<int> bn = pcoh.persistent_betti_numbers(0, 0);
        copy(bn.begin(), bn.end(), outp);
    }
}
