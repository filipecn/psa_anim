<div align="center">

# Digital Animation of Powder-Snow Avalanches

<h2>ACM Transactions on Graphics - SIGGRAPH 2025</h2>

[Filipe Nascimento](https://filipecn.dev/), [Fabricio S. Sousa](https://sites.google.com/icmc.usp.br/fssousa/home), [Afonso Paiva](https://sites.google.com/icmc.usp.br/apneto/)

University of Sao Paulo - (ICMC-USP)


### [Projectpage](https://filipecn.github.io/psa_anim/) · [Paper](https://filipecn.github.io/psa_anim/static/pdf/psa_anim.pdf) · [Video](https://www.youtube.com/watch?v=rHvtYA-lLIk&feature=youtu.be)

</div>

## Introduction

TODO

## C++ library Dependencies

The core code provides many tools written in C++ that are utilized by the 
several `bash`/`python` scripts in order to simulate and produce the final 
animation data. These tools rely on some external libraries, which are 
handled automatically by `cmake`, except for 

[OpenVDB](https://github.com/AcademySoftwareFoundation/openvdb)

[Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)

which are assumed to be available in the system before compilation. 
In case you have troubles with `OpenVDB` or `Eigen` during cmake/compilation,
you can try to modify the way they are included in the main `CMakeLists.txt`:
```
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenVDB REQUIRED)
```

## OpenFOAM Dependency

[OpenFOAM](https://www.openfoam.com/) is the framework we use to implement
our model. We particularly use the version 
[v2206](https://www.openfoam.com/news/main-news/openfoam-v2206). 
But there are newer versions that may work as well.

## Build

The compilation process utilizes the usual `cmake`/`make` ritual.

```
git clone https://github.com/filipecn/psa_anim.git
cd psa_anim 
mkdir build
cd build
cmake ..
make -j20
```

You will also need to compile the `OpenFOAM` solvers through `OpenFOAM`'s compiler `wmake`.

```
# faSavageHutterFoam
cd ext/faSavageHutterFoam
wmake

# pslFoam
cd pslFoam
wmake
```

## Usage

With both `OpenFOAM` solvers and the c++ tools (`tools` folder), you may utilize the 
many scripts found inside `tutorials` and `scripts` folders. In order to make these 
scripts visible in you shell, run 
```
source /path/to/psa_anim/build/psa_anim_variables.sh
```

For executing the full workflow, 
the most convenient way is to setup your case files like the tutorial 
cases and execute the `run_tutorial.sh` script. For example:

```
$PSA_ANIM_TUTORIALS/run_tutorial.sh --name niobe
```

This script accepts many arguments for configuring the simulation, such as duration, 
simulation parameters, output folders, etc. 
Please check the documentation for more details.

## Citation

TODO
