# Compiling CaNS

For most systems, CaNS can be compiled from the root directory with the following commands `make libs && make`, which will compile the 2DECOMP&FFT/cuDecomp libraries and CaNS. `make clean` clears the CaNS build files, `make libsclean` clears the 2DECOMP/cuDecomp builds, and `make allclean` clears both.

The `Makefile` in root directory is used to compiled the code, and is expected to work out-of-the-box for most systems. The `build.conf` file in the root directory can be used to choose the Fortran compiler (MPI wrapper), a few pre-defined profiles depending on the nature of the run (e.g., production vs debugging), and pre-processing options:

```shell
#
# compiler and compiling profile
#
FCOMP=GNU          # options: GNU, NVIDIA, INTEL
FFLAGS_OPT=1       # for production runs
FFLAGS_OPT_MAX=0   # for production runs (more aggressive optimization)
FFLAGS_DEBUG=0     # for debugging
FFLAGS_DEBUG_MAX=0 # for thorough debugging
#
# defines
#
SINGLE_PRECISION=0 # perform the whole calculation in single precision
GPU=0              # GPU build
```

In this file, `FCOMP` can be one of `GNU` (`gfortran`), `INTEL` (`ifort`), `NVIDIA` (`nvfortran`), or `CRAY` (`ftn`); the predefined profiles for compiler options can be selected by choosing one of the `FFLAGS_*` option; finer control of the compiler flags may be achieved by building with, e.g., `make FFLAGS+=[OTHER_FLAGS]`, or by tweaking the profiles directly under `configs/flags.mk`. Similarly, the library paths (e.g., for *FFTW*) may need to be adapted in the `Makefile` (`LIBS` variable) or by building with `make LIBS+='-L[PATH_TO_LIB] -l[NAME_OF_LIB]'`. Finally, the following pre-processing options are available:

 * `SINGLE_PRECISION` : calculation will be carried out in single precision (the default precision is double)
 * `GPU`              : enable GPU accelerated runs (requires the `FCOMP=NVIDIA`)

Typing `make libs` will build the 2DECOMP&FFT/cuDecomp libraries; then typing `make` will compile the code and copy the executable `cans` to a `run/` folder; `make run` will also copy the default input files `*.in` under `src/` to the same `run/` folder.

Note that cuDecomp needs to be dynamically linked before performing a GPU run. To do this, one should update the `LD_LIBRARY_PATH` environment variable as follows (from the root directory):
```shell
export LD_LIBRARY_PATH=$PWD/dependencies/cuDecomp/build/lib:$LD_LIBRARY_PATH
```

Finally, the choice of compiler `FCOMP` (see `configs/flags.mk`), and profile flags `FFLAGS_*` (see `configs/flags.mk`) can easily be overloaded, for instance, as: `make FC=ftn FFLAGS=-O2`. Linking and Include options can be changed in `configs/libs.mk`. The default `build.conf` files and `*.mk` files are created from `configs/defaults` at the first compilation.
