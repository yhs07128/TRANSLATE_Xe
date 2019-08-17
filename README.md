# LArCADe
Simulation of electron transport in LAr under non-uniform electric fields for LArCADe at FNAL.

Basis of simulation outlined in Wojcik, M., & Tachiya, M. (2002). Electron transport and electron–ion recombination in liquid argon: Simulation based on the Cohen–Lekner theory. Chemical Physics Letters, 363(3-4), 381-388. doi:10.1016/s0009-2614(02)01177-6

## Getting Started
Running the executable for LArCADe will walk you through the parameters to be specified before the simulation begins. The resulting data is exported to the py/simulation-runs folder. Examples of studies carried out with the produced simulations can be found at 

Gas cross-section data is from the BSR database, www.lxcat.net, retrieved on June 11, 2019 and Rejoub, R., Lindsay, B. G., & Stebbings, R. F. (2002). Determination of the absolute partial and total cross sections for electron-impact ionization of the rare gases. Physical Review A, 65(4). doi:10.1103/physreva.65.042713. Liquid cross-section data comes from the Wojcik, M., & Tachiya, M. paper mentioned at the top.

### Prerequisites
You'll need Numpy, Scipy, and Matplotlib. Jupyter is required to view the preprocessing notebooks. CMake makes for an easy install.

### Installing
Navigate to the project directory in terminal, then type

```
mkdir build
cd build
cmake ..
make
```