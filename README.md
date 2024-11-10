# Diffusive MPPIC (DMPPIC) Project

## What the Project Does

MC is a new communication paradigm that uses molecules/particles as the information carrier. We aim to simulate the channel impulse response (CIR) using the MPPICFoam solver in OpenFOAM. However, we found that the particles in MPPICFoam do not diffuse in the absence of external flow.

This project enhances the MPPICFoam solver in OpenFOAM by introducing diffusion effects. The improved version is called **Diffusive MPPIC (DMPPIC)**. In DMPPIC, we implement the diffusion effect of particles by adding random walk calculations into the source code.

**Most of the code used in this project is derived from ESI-OpenFOAM**

https://develop.openfoam.com/Development/openfoam

## How Users Can Get Started with the Project

### 1. Setup OpenFOAM-v2312 Environment

#### 1.1 Download OpenFOAM

If you are using Ubuntu, follow these steps:

```bash
sudo apt update -y
sudo apt upgrade -y
sudo apt install build-essential make gdb cmake -y
curl https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get install openfoam2312-default
sudo apt install openfoam-selector
```

### 1.2 Setup Environmental Variables
Open the .bashrc file in a text editor:
``` bash
nano ~/.bashrc
```
Add the following lines to the end of the file:
```bash
export WM_PROJECT_DIR="/usr/lib/openfoam/openfoam2312"
export FOAM_APPBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/bin"
export FOAM_LIBBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/lib"
source $WM_PROJECT_DIR/etc/bashrc
```
Save and close the file, then source the updated .bashrc:
``` bash
source ~/.bashrc
```
Verify the installation by running a test case:
```bash
foamInstallationTest -full incompressible/simpleFoam/pitzDaily
```

### 1.3 Compile Additional Library
``` bash
cd $WM_PROJECT_DIR/applications/solvers/lagrangian/DPMFoam/
sudo bash -c "source $WM_PROJECT_DIR/etc/bashrc && ./Allwmake"
```

### 1.4 Copy the Library Files to Project
``` bash
cp $WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/lib/* -r ./lib/
cp $WM_PROJECT_DIR/applications/solvers/lagrangian/DPMFoam/DPMTurbulenceModels -r ./lib/
```

### 1.5 Build Project
``` bash
cd /path/to/your/project
mkdir build && cd build
cmake ..
make
```

## 2. Run Simulation Using Compiled DMPPICFoam Solver
After successfully building the project, you can run simulations using the compiled DMPPICFoam solver. Refer to the Usage Guide for detailed instructions on setting up and executing simulations.
