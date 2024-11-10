# Diffusive MPPIC (DMPPIC) Project

## What the Project Does

MC is a new communication paradigm that uses molecules/particles as the information carrier. We aim to simulate the channel impulse response (CIR) using the MPPICFoam solver in OpenFOAM. However, we found that the particles in MPPICFoam do not diffuse in the absence of external flow.

This project enhances the MPPICFoam solver in OpenFOAM by introducing diffusion effects. The improved version is called **Diffusive MPPIC (DMPPIC)**. In DMPPIC, we implement the diffusion effect of particles by adding random walk calculations into the source code.

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

#### 1. Open the .bashrc file in a text editor:
nano ~/.bashrc

#### Add the following lines to the end of the file:
```bash
export WM_PROJECT_DIR="/usr/lib/openfoam/openfoam2312"
export FOAM_APPBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/bin"
export FOAM_LIBBIN="$WM_PROJECT_DIR/platforms/linux64GccDPInt32Opt/lib"
source $WM_PROJECT_DIR/etc/bashrc
```


