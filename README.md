# LArCADe
Simulation of electron transport in LAr for LArCADe at FNAL.

Basis of simulation outlined in Wojcik, M., & Tachiya, M. (2002). Electron transport and electron–ion recombination in liquid argon: Simulation based on the Cohen–Lekner theory. Chemical Physics Letters, 363(3-4), 381-388. doi:10.1016/s0009-2614(02)01177-6

## Getting Started
Running the executable for LArCADe will walk you through the parameters that need to be specified before simulation begins. To view the result, type

```
python grapher.py
```

while still in the build folder. If you specify

```
python grapher.py -m
```

the grapher will plot the electron mobility instead of its drift velocity. You can plot individual groups of simulation runs with

```
python grapher.py -g "BEGINNING OF FILE NAME"
```

For example,

```
python grapher.py -g "1D - 2000V"
```

Cross-sections can be plotted with by navigating to the py folder and typing

```
python cross_sections.py -d
```

Other arguments are -s (to save the cross-section data in text files for use in LUTs) and -p (to plot the cross-sections multiplied by their indexing velocities).

Gas cross-section data is from the BSR database, www.lxcat.net, retrieved on June 11, 2019 and Rejoub, R., Lindsay, B. G., & Stebbings, R. F. (2002). Determination of the absolute partial and total cross sections for electron-impact ionization of the rare gases. Physical Review A, 65(4). doi:10.1103/physreva.65.042713. Liquid cross-section data comes from the Wojcik, M., & Tachiya, M. paper mentioned at the top.

### Prerequisites
You'll need Numpy, Scipy, and Matplotlib. CMake makes for an easy install.

### Installing
The usual. Navigate to the project directory in terminal, then type

```
cd build
cmake ..
make
```