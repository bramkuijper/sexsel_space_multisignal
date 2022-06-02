#ifndef _SIMULATION_HPP_
#define _SIMULATION_HPP_

#include <vector>
#include <unordered_set>
#include <random>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include "patch.hpp"
#include "parameters.hpp"


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

        int mean_survivors_f = 0;
        int mean_survivors_m = 0;

    public:

        Parameters parms;

        std::vector<Patch> metapop;

        Simulation(Parameters const &parameters);

        // function that initializes all the patches
        // and their contents (breeders, environments)
        void initialize_patches();

        // initialize the data files
        void initialize_output_file();

        // function governing survival of carrying capacity
        void survive_carrying_capacity();

        void male_male_competition();

        void female_choice();

        void offspring_production();

        double carrying_capacity(double const environment_location);

        // sample k individuals out of n
        void sample_k_out_of_n(
                int const N
                ,int const k
                ,std::unordered_set<int> &individuals_sampled);


        void write_data();
        void write_parameters();
}; // end Simulation class definition

#endif
