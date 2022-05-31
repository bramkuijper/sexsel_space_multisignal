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

    // initial value of the local adaptation locus
    double v_env_init = 0.5;

    // number of males and females
    // per patch
    unsigned nm = 10;
    unsigned nf = 10;

    // number of patches in the population
    unsigned Npatches = 50;

    // dispersal
    double d = 0.5;


    // baseline carrying capacity
    double b = 0.0;

    // scalar to set the maximum carrying capacity
    // see M'Gonigle's supp
    double k0 = 1.0;


    // maximum time 
    long unsigned max_time = 100;

    // spatial variance in carrying capacity
    double sigma_k = 0.5;

    // the prefix of the file name
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
    } // end currentDateTime

    // constructor function of the parameter object
    // I do this to initialize the output 
    // file to something sensible
    Parameters() :
        base_name{"sim_sexsel_space_" + currentDateTime()} // upon construction, initialize the file name
    {}
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
        long unsigned time_step;
    public:

        Parameters parms;

        std::vector<Patch> metapop;

        Simulation(Parameters const &parameters);

        // function that initializes all the patches
        // and their contents (breeders, environments)
        void initialize_patches();

        // initialize the data files
        void initialize_output_file();

        // function governing offspring production and survival
        void offspring_production_and_survival();


        void male_male_competition();

        void female_choice();

        double carrying_capacity(double const environment_location);

        void sample_k_from_range(
            int const N, 
            int const k, 
            std::vector <int> &sampled_vector)
{

}; // end Simulation class definition

#endif
