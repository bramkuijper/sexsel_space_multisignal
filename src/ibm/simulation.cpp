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


// initialize the data files used to output the simulation results
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
        //
        // see individual.hpp and individual.cpp for definitions
        // of the individual
        Individual a_individual(
                    parms.p_loc_init
                    ,parms.t_loc_init
                    ,parms.p_comp_init
                    ,parms.t_comp_init
                    ,parms.v_env_init
                );

        // make a patch and then initialize it
        // see patch.cpp and patch.hpp
        Patch current_patch(
                parms.nm // number males
                ,a_individual // the standard male with which to initialize all males
                ,parms.nf // number of females
                ,a_individual // the standard female with which to initialize all females
                ,uniform(rng_r)
                );

        // give patches an environmental variable
        // along a 0 to 1 uniform distribution
        current_patch.coordinate = uniform(rng_r);

        // now create males and females

        metapop.push_back(current_patch);



    } // end for patch_idx
} // end Simulation::initialize_patches()


// environmental survival of juveniles

// production of offspring by females and survival
// based on that
void Simulation::offspring_production_and_survival()
{
    // auxiliary variable containing the carrying capacity of a patch
    double k;

    int nmales_sampled;
    int nfemales_sampled;

    std::vector <int> males_sampled;
    std::vector <int> females_sampled;

    // go through all the patches and perform births
    for (unsigned patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        // survival probability based on carrying capacity
        k = carrying_capacity(metapop[patch_idx].coordinate);

        // make a sampler out of a binomial distribution
        std::binomial_distribution<int> male_sampler{
            metapop[patch_idx].males.size(), k};

        // sample number of males
        nmales_sampled = males_sampler(rng_r);

    } // end for unsigned patch_idx
} // Simulation::offspring_production_and_survival

// Floyd's algorithm to sample k individuals out of a 
// list of N 
void Simulation::sample_k_from_range(
        int const N, 
        int const k, 
        std::vector <int> &sampled_vector)
{
    sampled_vector.clear();

    for (int r = N - k; r < N; ++r)
    {
        int v = std::uniform_int_distribution<int>(0,r)(rng_r);

        if (!elems.insert(v))
        {


    } // end for


}


// calculate unidimensional carrying capacity according to eq S1
// from MGonigle et al 2012 Nature 484: 506
double Simulation::carrying_capacity(double const coordinate)
{
    return(parms.k0 * (
                parms.b + exp(-(coordinate - 0.5)*
                    (coordinate - 0.5)/(2 * parms.sigma_k * parms.sigma_k))));
}

// competition among males dependent on hawk dove game
void Simulation::male_male_competition()
{}

void Simulation::female_choice()
{}
