## Sexual selection, varying densities and different signal types

### Reading the code
If you want to start reading the code, best to start with `simulation.cpp` and the corresponding header file `simulation.hpp`. This contains the main part of the code, as well as various functions that underlie the life cycle. The file `simulation.hpp` also contains a listing of all the relevant parameters used in the simulation. 

### Running and compiling the simulation
In order to run the thing, it is assumed you run a unix-like terminal (e.g., [msys2](https://www.msys2.org/) on Windows, or Terminal on mac os). You need to have installed `cmake` and a compiler.

In the terminal, run the following commands:
```
git clone https://github.com/bramkuijper/sexsel_space_multisignal.git 
cd sexsel_space_multisignal/src/ibm
cmake -S . -B build/
cmake --build build
cd build
./simulation_main
```

### Analysing the simulation
After the simulation is finished, you should see an output file named `sim_sexsel_space_YYYY_mm_dd`. If one opens the file using a text editor (like notepad++ or textmate), you see a collection of comma-separated values. This file contains all the statistics produced over the course of the simulation.

Using the `R`-package `bramkuijper:simulation.utils` one can start generating summary files of multiple simulations. 
