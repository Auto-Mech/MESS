[![C/C++ CI](https://github.com/keceli/MESS/actions/workflows/build.yml/badge.svg)](https://github.com/keceli/MESS/actions/workflows/build.yml)

# MESS

Master Equation System Solver

The primary purpose of the program is to calculate temperature and pressure dependent
rate coefficients for complex-forming reactions via solution of the one-dimensional master equation. Ancillary calculations of various quantities (e.g., stabilization probabilities for microcanonical initial distributions, microcanonical rate constants, partition functions, and related thermochemical information, time dependent propagation of species populations, etc.) are also available. 

### Requirements ###

To build MESS, the following libraries are required:  
BLAS  
LAPACK  
SLATEC <https://github.com/Auto-Mech/SLATEC>  
MPACK <https://github.com/Auto-Mech/MPACK>


### Direct Installation using Conda

The most direct way to install the code is through the conda package manager.
If you have conda installed,  
(1) activate an environment in you wish to use to install MESS, and  
(2) run the install command:
```
conda install -c auto-mech mess
```

If you do not have a preferred Conda environment set up, an empty environment with no packages can be created and activated with the following commands
```
conda create --name myenv
conda activate myenv
```
where `myenv` should be replaced with your preferred name for the environment.

Alternatively, we also recommend building our own pre-set Auto-Mech environment, which includes MESS and all its dependencies. This environment can be created and activated with the commands:
```
conda env create auto-mech/amech-env
conda activate amech-env
```

If your Conda commands are not functioning, you may need to iniliatize Conda via the command
```
. /path/to/conda.sh
```
which places Conda executables in your PATH. The specific location of conda.sh depends on the Conda install.

If you do not have Conda, it can be installed using the shell script
`debug/install-conda.sh`.


### Building from source using Conda environment for dependencies

To build the code from source for development or debugging purposes, first
create a conda environment with the necessary dependencies as follows:
```
conda env create -f environment.yml
```
which will create the `mess-env` environment.
You can then activate the environment and build the code as follows:
```
conda activate mess-env
bash debug/build.sh
```
To put the MESS executables in your PATH, you can then run
```
. debug/fake-install.sh
```
Note that the above command does not **permanently** alter your PATH, it only affects PATH for the current login session.


### Building from source without Conda

This is not the advised way to install, since the user will have to deal with their specific system setup.

Download SLATEC and MPACK from GitHub, and install them in a location that your system can find them. 

With SLATEC and MPACK installed, run build.sh, which uses cmake to compile MESS.
```
bash build.sh
```

Note that the results of the `make install` command in build.sh will depend on your system setup.

### Building from Source (CMake)

This project uses CMake for building from source. Ensure you have CMake (version 3.16 or higher), a C++ compiler (supporting C++11), and a Fortran compiler installed.

1.  **Clone the repository (if you haven't already):**
    ```bash
    git clone https://github.com/keceli/MESS.git
    cd MESS
    ```

2.  **Create a build directory and navigate into it:**
    ```bash
    mkdir build
    cd build
    ```

3.  **Configure the project using CMake:**
    Run CMake from the `build` directory, pointing to the parent directory (where the main `CMakeLists.txt` is located).
    ```bash
    cmake ..
    ```
    This command prepares the build system. By default, it will try to find system libraries and download/build dependencies if they are not found.

    **Build Options (passed to the `cmake ..` command):

    *   `AUTO_DOWNLOAD_DEPENDENCIES` (Default: `ON`): Automatically download and build missing dependencies. To disable, use:
        ```bash
        cmake -DAUTO_DOWNLOAD_DEPENDENCIES=OFF ..
        ```
    *   `USE_SYSTEM_LIBS` (Default: `ON`): Try to use system libraries before downloading. To disable (and force download/build of all dependencies), use:
        ```bash
        cmake -DUSE_SYSTEM_LIBS=OFF ..
        ```
    *   You can combine options: `cmake -DAUTO_DOWNLOAD_DEPENDENCIES=OFF -DUSE_SYSTEM_LIBS=OFF ..`

4.  **Compile the project:**
    After CMake configuration is successful, run `make` (or your chosen build tool like `ninja`) in the `build` directory.
    ```bash
    make -jN # Replace N with the number of parallel jobs you want to use, e.g., make -j4
    ```
    This will compile the MESS executables (mess, mess-v2, messpf, messabs, messsym) and place them in the `build` directory.

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
