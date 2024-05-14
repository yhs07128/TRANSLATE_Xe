Parameter tuned version of TRANSLATE(https://github.com/davidc1/TRANSLATE) to work for LXe.

## Getting Started
Running the executable for LArCADe will walk you through the parameters to be specified before the simulation begins. The resulting data is exported to the py/simulation-runs folder. 

### Prerequisites
You'll need Numpy, Scipy, and Matplotlib. Jupyter is required to view the preprocessing notebooks. CMake makes for an easy install.

### Installing
Open the Jupyter notebooks in the preprocessing directory to generate the lookup tables for cross-sections and non-uniform electric fields. After doing so, set the necessary parameters in the Constants.h file. Then, navigate to the project directory in terminal and type

```
mkdir build
cd build
cmake ..
make
```

### Setup
To run the code type `./larcade` and follow the prompts. To change running configuration edit the `srcs/Constants.h` file and specifically update the value of `n`, the density of argon atoms, `is_diff` and `is_gas` to choose the simulation mode.

Configurable constants:

```
// Liquid number density [cm^-3]                                                                                                                       
const double n = 1.35e22;       // Number density of argon atoms

// Editable parameters                                                                                                                                 
const bool is_diff = false;     // Use differential cross-sections?                                                                    
const bool is_gas = false;      // Is this in gas? Or liquid? Make sure to update the number density 
```

Default liquid mode (integrated xsec): `n = 1.35e22`, `is_diff = false`, and `is_gas = false`.

### Simulation Output
Information for each simulated electron is stored in text files. Each iteration of the simulation for which output is stored is saved as a new row, and information in each row is arranged as described below, separated by commas:

time (ns), time-step (ns), x coordinate (um), y coordainte (um), z coordainte (um), kinetic energy (eV), drift velocity (m/s), scattering angle (degrees), distance since last interaction (meters), interaction category, total ionizations initiated

