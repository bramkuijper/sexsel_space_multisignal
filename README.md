## Sexual selection, varying densities and different signal types

### Reading the code
If you want to start reading the code, best to start with `simulation.cpp` and the corresponding header file `simulation.hpp`. This contains the main part of the code, as well as various functions that underlie the life cycle. The file `simulation.hpp` also contains a listing of all the relevant parameters used in the simulation. 

### Running and compiling the simulation
In order to run the thing, it is assumed you run a unix-like terminal (e.g., [msys2](https://www.msys2.org/) on Windows, or Terminal on mac os). You need to have installed `cmake`, `git` and a compiler like `clang` (mac os) or `gcc`.

In the terminal, run the following commands to download and build the simulation
```
# download the code to the local directory sexsel_space_multisignal
git clone https://github.com/bramkuijper/sexsel_space_multisignal.git 

# enter directory containing the source code
cd sexsel_space_multisignal/src/ibm

# let cmake know where the source code can be found and in which 
# directory to build the executable 
cmake -S . -B build/

# then build the thing
cmake --build build
```

After this is all done, you can run the simulation as follows
```
cd build
./simulation_main 
```

### Analysing the simulation
After the simulation is finished, you should see an output file named `sim_sexsel_space_YYYY_mm_dd`. If one opens the file using a text editor (like notepad++ or textmate), one finds a collection of comma-separated values. This file contains all the statistics produced over the course of the simulation.

Using the `R`-package `bramkuijper:simulation.utils` one can start generating summary files of multiple simulations. 


