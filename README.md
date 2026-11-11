[![Build and test](https://github.com/Auto-Mech/MESS/actions/workflows/build.yml/badge.svg)](https://github.com/Auto-Mech/MESS/actions/workflows/build.yml)

# MESS

Master Equation System Solver

The primary purpose of the program is to calculate temperature and pressure dependent
rate coefficients for complex-forming reactions via solution of the one-dimensional master equation. Ancillary calculations of various quantities (e.g., stabilization probabilities for microcanonical initial distributions, microcanonical rate constants, partition functions, and related thermochemical information, time dependent propagation of species populations, etc.) are also available. 

## Installation

The [`automech`](https://prefix.dev/channels/automech) channel on prefix.dev provides two conda packages with the same executables (`mess`, `messpf`, `messabs`).
Install only one of them.

`mess-static` has fully static executables that run on any x86-64 Linux, with no runtime dependencies:
```
conda install -c https://repo.prefix.dev/automech mess-static
```
In a [Pixi](https://pixi.sh) workspace:
```
pixi workspace channel add https://repo.prefix.dev/automech
pixi add mess-static
```
`mess` is linked against the conda-forge libraries (OpenBLAS, GSL), so it also needs the `conda-forge` channel:
```
conda install -c https://repo.prefix.dev/automech -c conda-forge mess
```
The static executables are also attached to each [GitHub release](https://github.com/Auto-Mech/MESS/releases) as a tarball, for systems without conda.
The changes in each release are listed in [`CHANGELOG.md`](CHANGELOG.md).

The `auto-mech` channel on anaconda.org keeps the `mess-static` packages up to version 2026.2.4 and gets no new ones.

## Building from source

The build uses [Pixi](https://pixi.sh), with all dependencies from conda-forge:
- ILP64 BLAS/LAPACK: `openblas-ilp64`
- GSL
- OpenMP
- SLATEC: the routines MESS needs are vendored in `external/SLATEC`

From the repository root:
```
pixi run test             # build MESS (linked against conda-forge libraries) and run the examples
pixi run -e static test   # the same for the fully static executables
```
The executables (`mess`, `messpf`, `messabs`) are in `.pixi/envs/default/bin/` and `.pixi/envs/static/bin/`.
`pixi shell` puts the default ones on your `PATH`.

### Plain CMake

`CMakeLists.txt` is a standard CMake project. It needs C, C++ and Fortran compilers, an ILP64 (64-bit integer) BLAS/LAPACK, GSL and OpenMP.
```
cmake -B build -DBLAS_LIBRARIES=/path/to/libopenblas64_.so -DLAPACK_LIBRARIES=/path/to/libopenblas64_.so -DMESS_BLAS_SUFFIX=64_
cmake --build build
```
Options:
- `BLAS_LIBRARIES`, `LAPACK_LIBRARIES`: the ILP64 library. These can be left out if CMake's FindBLAS finds it (e.g. `-DBLA_VENDOR=Intel10_64ilp` for MKL).
- `MESS_BLAS_SUFFIX`: the suffix on the library's BLAS/LAPACK symbol names, if any. For example, `64_` for OpenBLAS built with `SYMBOLSUFFIX=64_`, like conda-forge's `libopenblas64_.so`. Check with `nm -D <library> | grep dsyev`. The static archive `libopenblas64_.a` has no suffix.
- `MESS_STATIC` (default `OFF`): link fully static executables.

## Releasing

`mess` and `mess-static` share one version, `YYYY.M.N`, where `N` counts the releases within the month, starting from 0.
To release, from an up-to-date `master` of this repository:
1. List the changes to the source since the last release:
   ```
   git log --no-merges $(git describe --tags --abbrev=0)..HEAD -- src examples
   ```
2. Describe them under `[Unreleased]` in `CHANGELOG.md`, and commit.
3. Run `pixi run release YYYY.M.N` (for example, `pixi run release 2026.11.0`).
   This updates the version and the changelog, then commits, tags and pushes.

The tag starts the release workflow (`.github/workflows/release.yml`), which builds and tests both packages, publishes them to the `automech` channel, and makes a GitHub release.
`pixi run version` shows the current version.

## Reference

See Y. Georgievskii, J. A. Miller, M. P. Burke, and S. J. Klippenstein,
Reformulation and Solution of the Master Equation for Multiple-Well Chemical
Reactions, J. Phys. Chem. A, 117, 12146-12154 (2013).

## Acknowledgment

This work was supported by the U.S. Department of Energy, Office of Basic Energy
Sciences, Division of Chemical Sciences, Geosciences, and Biosciences under DOE
Contract Number DE-AC02-06CH11357 as well as the Exascale Computing Project
(ECP), Project Number: 17-SC-20-SC.  The ECP is a collaborative effort of two
DOE organizations, the Office of Science and the National Nuclear Security
Administration, responsible for the planning and preparation of a capable
exascale ecosystem including software, applications, hardware, advanced system
engineering, and early test bed platforms to support the nation's exascale
computing imperative. 

## Notice

Copyright (c) 2018 Yuri Georgievski (ygeorgi@anl.gov), Stephen J.
Klippenstein (sjk@anl.gov), and Argonne National Laboratory.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
