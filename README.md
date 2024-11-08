What the project does
    MC is a new communication paradim that uses molecules/particles as the information carrier.
    We want to simulation the channel impulse response (CIR) using the MPPICFoam in OpenFOAM.
    However, we found the particles in MPPICFoam do not diffuse if there is not external flow.
    This project improves the MPPICFoam solver in OpenFOAM slightly, we called the improved version as Diffusive MPPIC (DMPPIC).
    In DMPPIC, we implement diffusion effect of particles by adding the radom walk calculations into the source code.

How users can get started with the project
1. Setup OpenFOAM-v2312 environment.

1.1 Download OpenFOAM
If use Ubuntu:
$ sudo apt update -y
$ sudo apt upgrade -y
$ sudo apt install build-essential make gdb cmake -y
$ curl https://dl.openfoam.com/add-debian-repo.sh | sudo bash
$ sudo apt-get install openfoam2312-default
$ sudo apt install openfoam-selector

1.2 Setup environmental variables
$ nano ~/.bashrc
Add content below:
export WM_PROJECT_DIR="/usr/lib/openfoam/openfoam2312"
export FOAM_APPBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/bin"
export FOAM_LIBBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/lib"
source $WM_PROJECT_DIR/etc/bashrc
$ source ~/.bashrc
$ foamInstallationTest -full incompressible/simpleFoam/pitzDaily

1.3 Compile additional library
$ cd $WM_PROJECT_DIR/applications/solvers/lagrangian/DPMFoam/
$ sudo bash -c "source $WM_PROJECT_DIR/etc/bashrc && ./Allwmake"

1.3 Copy the library files to project
cp $WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/lib/* -r ./lib/
cp $WM_PROJECT_DIR/applications/solvers/lagrangian/DPMFoam/DPMTurbulenceModels -r ./lib/

1.4 Build project
cd to the project directory
$ mkdir build && cd build
$ cmake ..
$ make

2. Run simulation using compiled DMPPICFoam solver