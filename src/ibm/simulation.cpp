#include <unordered_set>
#include "simulation.hpp"
#include "parameters.hpp"
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
        survive_carrying_capacity();

        male_male_competition();

        // female choice based on
        // 1. indirect benefits of choice where
        // genetic quality is the result of local
        // adaptation 
        // 2. male-male competition
        female_choice();

        offspring_production();


        if (time_step % parms.output_interval == 0)
        {
            write_data();
        }
    }

    write_parameters();
} // Simulation::Simulation()


void Simulation::write_parameters()
{
    output_file << std::endl
        << std::endl
        << "nm;" << parms.nm << std::endl
        << "nf;" << parms.nf << std::endl;
} // write_parameters()

void Simulation::write_data()
{
    output_file << time_step << ";"
        << mean_survivors_f << ";"
        << mean_survivors_m << ";";

    output_file << std::endl;
} // end Simulation::write_data()


// initialize the data files used to output the simulation results
void Simulation::initialize_output_file()
{
    output_file << "generation;mean_survivors_f;mean_survivors_m;mean_p_loc;mean_t_loc;mean_p_comp;mean_t_comp;mean_envt;mean_v_env;" << std::endl;
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

        metapop.push_back(current_patch);
    } // end for patch_idx
} // end Simulation::initialize_patches()


// environmental survival of juveniles

// survival based on carrying capacity
void Simulation::survive_carrying_capacity()
{
    // auxiliary variable containing the carrying capacity of a patch
    double k;

    int nmales_sampled;
    int nfemales_sampled;

    int mean_survivors_f = 0;
    int mean_survivors_m = 0;

    std::unordered_set<int> males_sampled;
    std::unordered_set<int> females_sampled;

    std::vector<Individual> surviving_males;
    std::vector<Individual> surviving_females;

    // go through all the patches and perform births
    for (unsigned patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        surviving_males.clear();
        surviving_females.clear();


        // survival probability based on carrying capacity
        k = carrying_capacity(metapop[patch_idx].coordinate);

        if (metapop[patch_idx].males.size() > 0)
        {
            // sample surviving males and females
            // how many males will survive?
            // make a binomial distribution 
            // with parameters 
            // 1) the number of males
            // 2) the probability of surviving, k 
            std::binomial_distribution<int> nmale_sample_dist(
                    metapop[patch_idx].males.size()
                    ,k);

            nmales_sampled = nmale_sample_dist(rng_r);

            sample_k_out_of_n(metapop[patch_idx].males.size() - 1
                    ,nmales_sampled
                    ,males_sampled);

            // loop through all males and have them survive
            for (std::unordered_set<int>::iterator males_iter = males_sampled.begin();
                    males_iter != males_sampled.end();
                    ++males_iter)
            {
                surviving_males.push_back(metapop[patch_idx].males[*males_iter]);
            }
        }

        if (metapop[patch_idx].females.size() > 0)
        {
            // now the females
            std::binomial_distribution<int> nfemale_sample_dist(
                    metapop[patch_idx].females.size()
                    ,k);

            nfemales_sampled = nfemale_sample_dist(rng_r);

            sample_k_out_of_n(metapop[patch_idx].males.size() - 1
                    ,nfemales_sampled
                    ,females_sampled);

            // loop through all males and have them survive
            for (std::unordered_set<int>::iterator females_iter = females_sampled.begin();
                    females_iter != females_sampled.end();
                    ++females_iter)
            {
                surviving_females.push_back(metapop[patch_idx].females[*females_iter]);
            }
        }

        // keep survivors
        metapop[patch_idx].males = surviving_males;
        metapop[patch_idx].females = surviving_females;

        mean_survivors_f += metapop[patch_idx].females.size();
        mean_survivors_m += metapop[patch_idx].males.size();
    } // end for unsigned patch_idx
      //
    mean_survivors_f /= metapop.size();
    mean_survivors_m /= metapop.size();
} // Simulation::survive_carrying_capacity

// calculate unidimensional carrying capacity according to eq S1
// from MGonigle et al 2012 Nature 484: 506
double Simulation::carrying_capacity(double const coordinate)
{
    return(parms.k0 * (
                parms.b + exp(-(coordinate - 0.5)*
                    (coordinate - 0.5)/(2 * parms.sigma_k * parms.sigma_k))));
}

// sample k out of n individuals
// this is floyd's algorithm, see
// https://stackoverflow.com/questions/28287138/c-randomly-sample-k-numbers-from-range-0n-1-n-k-without-replacement 
void Simulation::sample_k_out_of_n(
        int const N
        ,int const k
        ,std::unordered_set<int> &individuals_sampled)
{
    individuals_sampled.clear();

    for (int r = N - k; r < N; ++r)
    {
        int v = std::uniform_int_distribution<int>(0,r)(rng_r);

        // second tells you whether v was successfully inserted
        // into the unordered set elems. However, if v was already
        // in elems previously, it will not be inserted again (coz, sets)
        // and then second is false
        if (!individuals_sampled.insert(v).second)
        {
            individuals_sampled.insert(r);
        }
    }
}

// competition among males dependent on hawk dove game
void Simulation::male_male_competition()
{}

void Simulation::female_choice()
{
    // go through all the patches and perform choice
    for (std::vector<Patch>::iterator patch_iter = metapop.begin(); 
            patch_iter != metapop.end(); 
            ++patch_iter)
    {
        std::cout << patch_iter->females.size() << " ";
        std::cout << patch_iter->males.size() << " " << std::endl;
        for (std::vector<Individual>::iterator female_iter = patch_iter->females.begin(); 
                female_iter != patch_iter->females.end(); 
                ++female_iter)
        {
        }
    }
}

// produce a bunch of offspring
void Simulation::offspring_production()
{

}
