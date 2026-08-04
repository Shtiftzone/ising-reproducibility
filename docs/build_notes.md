# Build notes

## GUDHI dependency

The topological part of the code uses GUDHI in betti3d.cc.

Relevant includes:

- gudhi/Bitmap_cubical_complex.h
- gudhi/Persistent_cohomology.h

The current Makefile expects GUDHI headers to be available at:

src/simulation/gudhi/include/

For the current HPC setup, this can be provided by a symbolic link:

cd src/simulation
ln -s /data/home/mpx641/Ising/IsingModels/gudhi gudhi

The code also requires Boost, because GUDHI includes headers such as:

boost/intrusive/set.hpp

Current build status:

- GUDHI is found through the local symlink.
- Compilation currently fails because Boost headers are not found.
- The next step is to identify the required Boost module on Apocrita.

Before public release, the exact GUDHI and Boost versions should be identified and the installation instructions should be made independent of the local HPC path.

The CELL variants ising2d4-* and isingtr6-* are not built by default because they were not used in the analysis.
## Current working build

The simulation code builds successfully on Apocrita with:

module load boost/1.85.0-gcc-12.2.0
cd src/simulation
ln -s /data/home/mpx641/Ising/IsingModels/gudhi gudhi
make

The default Makefile builds only the executables used in the analysis. The CELL variants ising2d4-* and isingtr6-* are not built by default because they were not used in the analysis.

## Multihistogram build

The multihistogram analysis code is located in:

src/multihistogram/

It builds successfully on Apocrita with:

cd src/multihistogram
make

The resulting executables are build artifacts and are not tracked by Git.
