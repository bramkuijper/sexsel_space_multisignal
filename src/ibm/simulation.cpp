#include "simulation.hpp"
#include "patch.hpp"

//! construct a simulation object
//! @param  parameters  A parameter struct (see simulation.hpp)
//!                     that contains a collection of all the different
//!                     parameters that are needed to run the simulation
Simulation::Simulation(Parameters const &parameters) :
    rng_r(rd()) // initialize the random number generator with a random_device
    ,uniform(0.0,0.01) // initialize a uniform (0,1) distribution
    ,output_file(parameters.base_name.c_str())
    ,parms{parameters} // initialize the parameter member variable
{
    // initialize the data output file that will collect the data
    initialize_output_file();

    // initialize all the patches and their environment
    initialize_patches();

    for (time_step = 0; time_step < parms.max_time; ++time_step)
    {
        // environmental survival of juveniles
        offspring_production_and_survival();

        male_male_competition();

        // female choice based on
        // 1. indirect benefits of choice where
        // genetic quality is the result of local
        // adaptation 
        // 2. male-male competition
        female_choice();
    }
} // Simulation::Simulation()


// production of offspring by females and survival
// based on that
void Simulation::offspring_production_and_survival()
{
    for (unsigned patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {

    }
} // Simulation::offspring_production_and_survival

// initialize the data files
void Simulation::initialize_output_file()
{
    output_file << "mean_p_loc;mean_t_loc;mean_p_comp;mean_t_comp;mean_envt;mean_v_env;";
}// end Simulation::initialize_output_files()

// initialize all their patches 
// and their carrying capacities
void Simulation::initialize_patches()
{
    // loop through all patches and set carrying capacity
    for (int patch_idx = 0; patch_idx < parms.Npatches; ++patch_idx)
    {
        // make a default individual 
        // and use this to initialize patches
        Individual a_individual(
                    parms.p_loc_init
                    ,parms.t_loc_init
                    ,parms.p_comp_init
                    ,parms.t_comp_init
                    ,parms.venv
                );
        // make a patch and then initialize it
        Patch current_patch(
                parms.nm
                ,
                );

        // give patches an environmental variable
        // along a 0 to 1 uniform distribution
        current_patch.env = uniform(rng_r);

        // now create males and females

        metapop.push_back(current_patch);



    } // end for patch_idx
} // end Simulation::initialize_patches()
