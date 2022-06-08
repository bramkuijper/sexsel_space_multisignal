#include "parameters.hpp"
#include "simulation.hpp"


int main(int argc, char **argv)
{
    // change this bit of code if you want to have more parameters
    // the array argv cointains the command line parameters
    //
    Parameters params;
    params.max_time = atoi(argv[1]);
    params.d = atof(argv[2]);
    params.base_name = argv[3];

    Simulation the_sim{params};
}
