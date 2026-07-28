# CuCTMC
A set of parallel tools to analyze Continuous-Time Markov Chains (CTMCs) derrived from the Chemical Master Equation (CME).

This repository was created to support the Martinson Applied Project (MAP) "Stochastic Modeling of Host Immune Response to Viral Infections" at the University of Pittsburgh, led by Dr. David Swigon and Dr. Gilles Clermont. 
I created this repository independently, in my free time. It is not an official Pitt or MAP repository, and it has not been reviewed or endorsed by the project leaders.

## Requirements

Before you install this software, ensure that your system has the required hardware and software. The requirements include:
- An NVIDIA GPU with [compute capability](https://developer.nvidia.com/cuda/gpus) 7.5 or higher
- The NVIDIA [CUDA Toolkit](https://developer.nvidia.com/cuda/toolkit) (version 13.0 or higher)
- An NVIDIA [driver](https://www.nvidia.com/en-us/drivers/) that is compatible with your GPU and CUDA Toolkit version

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



## License

[Apache 2.0](https://choosealicense.com/licenses/apache-2.0/)
