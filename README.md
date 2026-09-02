# CuCTMC
A set of parallel tools to analyze Continuous-Time Markov Chains (CTMCs) derrived from the Chemical Master Equation (CME).

This repository was created to support the Martinson Applied Project (MAP) "Stochastic Modeling of Host Immune Response to Viral Infections" at the University of Pittsburgh, led by Dr. David Swigon and Dr. Gilles Clermont. 
I created this repository independently, in my free time. It is not an official Pitt or MAP repository, and it has not been reviewed or endorsed by the project leaders.

## Requirements

Before you install this software, ensure that your system has the required hardware and software. The requirements include:
- An NVIDIA GPU with [compute capability](https://developer.nvidia.com/cuda/gpus) 7.5 or higher
- The NVIDIA [CUDA Toolkit](https://developer.nvidia.com/cuda/toolkit) (version 13.0 or higher)
- An NVIDIA [driver](https://www.nvidia.com/en-us/drivers/) that is compatible with your GPU and CUDA Toolkit version

### A note on operating systems

This software is not, and has no plans to be, optimized for the Windows family of operating systems.

If you choose to compile this software for Windows, know that WDDM adds significant overhead to every CUDA driver call. In its current state, this codebase makes many such calls per sweep, which could mean the difference between seconds and minutes in execution time. 

This, along with my general distaste for the included bloatware, spyware, and forced AI integration of Windows operating systems, is why the repo will be built and tested solely on Linux machines for the foreseeable future.

## Installation

This software uses [CMake](https://cmake.org/download/) version 3.24 or higher as a build system generator.

Clone or download the repository. From the repository directory, run: 

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=<Debug|Release> -DINSTALL_MODELS=<ON|OFF>
cmake --build build --parallel $(nproc)
cmake --install build --prefix </path/to/install>
```

If you set `INSTALL_MODELS` to `ON`, the installation includes the example models in the installed directory.

### Uninstall

From the repository directory, run: 

```bash
xargs rm < build/install_manifest.txt
```

## Usage

This software is designed to be used from the command line. From the install directory, run: 

```bash
./CuCTMC [OPTIONS] <reaction file.json>
```

### Reaction Files

Reaction networks are given in terms of reactions and reactants, specified in JSON format. Each file must specify three things. First, specify the name of each reactant. For example:

```json
"reactants" : [
  "H", "V", "I", "D"
],
```

Then, specify the initial conditions. For example:


```json
"conditions" : [
  1000, 3, 0, 0
],
```

Lastly, specify the reactions. Each reaction must include a name, a nonnegative rate, and a list of reactants and products, and their multiplicities. For example:


```json
"reactions" : [
  {
    "name" : "gamma_VH", 
    "rate" : 0.0007, 
    "reactants" : { "V" : 1, "H" : 1 }, 
    "products" : { "I" : 1 }
  },
  {
    "name" : "gamma_V", 
    "rate" : 37.5, 
    "reactants" : { "I" : 1 }, 
    "products" : {"V" : 1, "I" : 1 }
  },
  {
    "name" : "a_I", 
    "rate" : 1.5, 
    "reactants" : { "I" : 1 }, 
    "products" : { "D" : 1 }
  },
  {
    "name" : "a_V", 
    "rate" : 1.7, 
    "reactants" : { "V" : 1 }, 
    "products" : {}
  },
  {
    "name" : "b_HD", 
    "rate" : 0.004, 
    "reactants" : { "H" : 1, "D" : 1 }, 
    "products" : { "H" : 2 }
  }
]
```

See [models](/models) for more examples on constructing reaction networks.

### CLI Options
The following command line flags are available. For options that permit swept parameters, all sets of parameters must have the same number of values.

```bash
-w <n>, --warps <n>
-w 6
```

Use this flag to set the number of warps. `n` is the number of warps. One warp is 32 threads in hardware terms. Each warp runs 32 paths. `n` must be a positive integer. You cannot sweep this parameter.

```bash
-T <x>, --tMax <x>
-T 50
```

Use this flag to set the maximum time of each sample path. `x` is the maximum time. `x` must be a positive float. You cannot sweep this parameter.

```bash
--rate <name> <values...>
--rate GammaV 15
--rate GammaV $(seq 1.5 1.5 150)
```

Use this flag to set the rate of one reaction. `name` is the name of the reaction. `values` is one rate, or a list of rates. Use one rate to keep the rate fixed. Use a list of rates to sweep the rate over many simulations.

```bash
--reactants <name> <values...>
--reactants Gamma 1 0 1
--reactants Gamma $(printf '%i 0 1' $(seq 1 10))
```

Use this flag to set the reactants of one reaction. `name` is the name of the reaction. `values` is one set of reactant counts, with one count per species. Use one set to keep the reactants fixed. Use many sets to sweep the reactants over many simulations. 

```bash
--products <name> <values...>
--products aI 0 35 0 1
--products aI $(printf '0 %i 0 1' $(seq 1 50))
```

Use this flag to set the products of one reaction. `name` is the name of the reaction. `values` is one set of product counts, with one count per species. Use one set to keep the products fixed. Use many sets to sweep the products over many simulations. 

```bash
--ic <values...>
--ic 1000 3 0 0
--ic $(printf '1000 %i 0 0' $(seq 1 30))
```

Use this flag to set the initial condition of the system. `values` is one set of reactant counts, with one count per species. Use one set to keep the initial condition fixed. Use many sets to sweep the initial condition over many simulations. Each value must be a non-negative integer.

```bash
-o <path>, --output <path>
-o ./results
```

Use this flag to set the output directory. `path` is the directory for the output files. You must set this flag if you use the `-s, --save` flag.

```bash
-t <x>, --tGrid <x>
-t 0.5
```

Use this flag to set the time grid spacing for saved trajectories. `x` is the time between two saved points on a trajectory. `x` must be a non-negative float. You must set this flag if you use the `-s, --save` flag.

```bash
-s <n>, --save <n>
-s 10
```

Use this flag to set the number of paths to save. `n` is the number of paths for which the program saves a trajectory. `n` must be a non-negative integer. You must also set the `-t, --tGrid` flag and the `-o, --output` flag.

## License

[Apache 2.0](https://choosealicense.com/licenses/apache-2.0/)
