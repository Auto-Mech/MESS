# MESS

[![Build Status](https://travis-ci.com/PACChem/MESS.svg)](https://travis-ci.com/PACChem/MESS)

## Installation

### Linux
- sudo apt-get install gcc gfortran libopenblas-dev liblapack-dev
- cmake .
- sudo make install

### Mac

The AppleClang compiler does not work with OpenMP, so you need to install GCC.
By default, `g++` will refer to the AppleClang version, so you will need to
feed the appropriate compiler to `cmake`. The procedure will look something
like this:
- brew update
- brew install gcc libomp
- cmake . -DCMAKE_CXX_COMPILER=g++-8
- sudo make install


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
