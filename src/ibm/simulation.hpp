#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include "patch.hpp"

// parameters used in the simulation with some default values
struct Parameters
{
    // initial value of the preference locus
    // based on local adaptation
    double p_loc_init = 0.0;
    // initial value of the ornament locus
    // based on local adaptation
    double t_loc_init = 0.0;
    // initial value of the preference locus
    // based on male-male competition
    double p_comp_init = 0.0;
    // initial value of the ornament locus
    // based on male-male competition
    double t_comp_init = 0.0;

    int Npatches = 100;
    double d = 0.5;

    long unsigned max_time = 50000;

    std::string base_name = "";

    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
    // this is a bit complicated but 
    // see https://stackoverflow.com/a/10467633 
    std::string currentDateTime() {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);

        // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
        // for more information about date/time format
        strftime(buf, sizeof(buf), "%Y_%m_%d_%H_%M_%S", &tstruct);

        return buf;
    }

    // constructor function of the parameter object
    // I do this to initialize the output 
    // file to something sensible
    Parameters() :
        base_name{"sim_sexsel_space_" + currentDateTime()}
//        base_name{"sim_sexsel_space_" + "hoi"}
    {
    }
};

// the simulation class which is the overarching
// structure of the simulation
class Simulation
{
    private:
        std::ofstream output_file;

        std::random_device rd;
        std::mt19937 rng_r;
        std::uniform_real_distribution<double> uniform;
        long unsigned time_step
    public:

        Parameters parms;

        std::vector<Patch> metapop;

        Simulation(Parameters const &parameters);

        void initialize_patches();
        void initialize_output_file();
};

#endif
