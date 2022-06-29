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
#include "individual.hpp"
#include "parameters.hpp"


// the simulation class which is the overarching
// structure of the simulation
class Simulation
{
    private:
        std::ofstream output_file;
        long unsigned time_step;
        
        // random device which is used to generate
        // proper random seeds
        std::random_device rd;

        // store the random seed
        // we need to store this so that we can output the 
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;
        
        // random number generator
        std::mt19937 rng_r;
        
        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;


        double mean_survivors_f = 0;
        double mean_survivors_m = 0;

        // number of juveniles produced each generation
        int njuveniles = 0;

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

        // females choose males and make kids
        void female_choice();

        double carrying_capacity(double const environment_location);

        double female_survival_probability(double envt, Individual &female_i);
        double male_survival_probability(double envt, Individual &male_i);

        double calculate_attractiveness(
                Individual &female
                ,Individual &male);

        void write_data();
        void write_parameters();

        void dispersal_and_replacement();
}; // end Simulation class definition

#endif
