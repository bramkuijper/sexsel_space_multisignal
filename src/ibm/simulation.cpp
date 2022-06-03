#include <unordered_set>
#include <cassert>
#include "simulation.hpp"
#include "parameters.hpp"
#include "patch.hpp"

//! construct a simulation object
//! @param  parameters  A parameter struct (see simulation.hpp)
//!                     that contains a collection of all the different
//!                     parameters that are needed to run the simulation
Simulation::Simulation(Parameters const &parameters) :
    rng_r(rd()) // initialize the random number generator with a random_device
    ,uniform(0.0,1.0) // initialize a uniform (0,1) distribution
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

// write the parameters to a file
void Simulation::write_parameters()
{
    output_file << std::endl
        << std::endl
        << "nm;" << parms.nm << std::endl
        << "nf;" << parms.nf << std::endl;
        << "ct;" << parms.ct << std::endl
        << "cp;" << parms.cp << std::endl;
        << "npatches;" << parms.Npatches << std::endl;
        << "d;" << parms.d << std::endl;
        << "k0;" << parms.k0 << std::endl;
        << "sigma_k;" << parms.sigma_k << std::endl;
} // write_parameters()

// write the stats to a file
void Simulation::write_data()
{
    // set up variables to take stats
    double mean_t = 0;
    double var_t = 0;

    double mean_t_conddep = 0;
    double var_t_condep = 0;

    double mean_p = 0;
    double var_p = 0;

    double mean_v = 0;
    double var_v = 0;
    
    double mean_s = 0;
    double var_s = 0;

    double p,t,tpr,v,s;

    // go through all the patches and perform choice
    for (std::vector<Patch>::iterator patch_iter = metapop.begin(); 
            patch_iter != metapop.end(); 
            ++patch_iter)
    {
        for (std::vector<Individual>::iterator female_iter = patch_iter->females.begin(); 
                female_iter != patch_iter->females.end(); 
                ++female_iter)
        {
            p = female_iter->p_loc[0] + female_iter->p_loc[1];
            t = female_iter->t[0] + female_iter->t[1];
            tpr = female_iter->t_conddep[0] + female_iter->t_conddep[1];
            v = female_iter->v_env[0] + female_iter->v_env[1];
            s = t + tpr * v;

        }


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

// male survival probability
double Simulation::male_survival_probability(double envt, Individual &male_i)
{
    // express ornament locus
    double t = male_i.t_loc[0] + male_i.t_loc[1];

    // express local adaptation locus
    double v = male_i.v_env[0] + male_i.v_env[1];

    // ornament based survival is exp(-ct * t^2) as in Iwasa et al 1991
    // Evolution 45: 1431
    //
    // this is a revealing conddepicap - see p1437 2nd column, last paragraph
    //
    // viability based survival is based on local adaptation, where a
    // good match with the local environment raises survival prospects
    // according to a bell shaped function exp(-vstr * (envt - v)^2)

    double k= carrying_capacity(envt);
    
    return(
            exp(-parms.ct * t * t) * 
            exp(-parms.vstrength * (envt - v) * (envt - v)) * k
            );
} // void Simulation::male_survival_probability(Individual &male_i)


double Simulation::female_survival_probability(double envt, Individual &female_i)
{
    // express ornament locus
    double p = female_i.p_loc[0] + female_i.p_loc[1];

    // express local adaptation locus
    double v = female_i.v_env[0] + female_i.v_env[1];

    double k= carrying_capacity(envt);

    // preference based survival is exp(-cp * p^2) as in Iwasa et al 1991
    // Evolution 45: 1431
    //
    // viability based survival is based on local adaptation, where a
    // good match with the local environment raises survival prospects
    // according to a bell shaped function exp(-vstr * (envt - v)^2)
    return(
            exp(-parms.cp * p * p) * 
            exp(-parms.vstrength * (envt - v) * (envt - v)) * k
            );
}

// survival based on carrying capacity
void Simulation::survive_carrying_capacity()
{
    // reset the stats that measures survival rates
    mean_survivors_f = 0;
    mean_survivors_m = 0;


    // vectors containing those males and females that survive
    // these have to be cleared when iterating over patches
    // otherwise survivors in one patch carry over to the next patch
    std::vector<Individual> surviving_males;
    std::vector<Individual> surviving_females;
    
    // aux variable containing the local envt 
    double envt;

    // go through all the patches and perform births
    for (unsigned patch_idx = 0; patch_idx < metapop.size(); ++patch_idx)
    {
        envt = metapop[patch_idx].coordinate;

        // survivors have to be cleared when starting to 
        // do survival business in the new patch
        surviving_males.clear();
        surviving_females.clear();

        for (std::vector<Individual>::iterator male_iter = metapop[patch_idx].males.begin(); 
                male_iter != metapop[patch_idx].males.end(); 
                ++male_iter)
        {
            if (uniform(rng_r) < male_survival_probability(envt, *male_iter))
            {
                surviving_males.push_back(*male_iter);
            }
        }

        for (std::vector<Individual>::iterator female_iter = metapop[patch_idx].females.begin(); 
                female_iter != metapop[patch_idx].females.end(); 
                ++female_iter)
        {
            if (uniform(rng_r) < female_survival_probability(envt, *female_iter))
            {
                surviving_females.push_back(*female_iter);
            }
        }

        metapop[patch_idx].females = surviving_females;
        metapop[patch_idx].males = surviving_males;

        mean_survivors_f += metapop[patch_idx].females.size();
        mean_survivors_m += metapop[patch_idx].males.size();

    } // end for unsigned patch_idx
    
    
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
        for (std::vector<Individual>::iterator female_iter = patch_iter->females.begin(); 
                female_iter != patch_iter->females.end(); 
                ++female_iter)
        {
            std::vector<double> attractiveness_values;

            // loop through males and make distribution of all male attractiveness values
            for (std::vector<Individual>::iterator male_iter = patch_iter->males.begin(); 
                    male_iter != patch_iter->males.end(); 
                    ++male_iter)
            {
                // calculate attractiveness and add this to a vector
                attractiveness_values.push_back(
                        calculate_attractiveness(*female_iter, *male_iter, patch_iter->coordinate)
                        );
            }

            // set up a probability distribution of all the males and choose from them
            // males with larger attractiveness values are more likely to be chosen
            std::discrete_distribution<int> choose_male_dist(attractiveness_values.begin()
                    ,attractiveness_values.end());

            // actually
            male_chosen = choose_male_dist(rng_r);

        } // end for female
    } // end for patch
} // end Simulation::female_choice()

double Simulation::calculate_attractiveness(Individual &female
        Individual &male)
{
    // express ornament
    double p = female.p_loc[0] + female.p_loc[1];
    double t = male.t_loc[0] + male.t_loc[1];
    double tprime = male.t_loc_conddep[0] + male.t_loc_conddep[1];
    double v = male.v_env[0] + male.v_env[1];

    // express male ornament as a condition-dependent trait
    // hence we assume s = t + tcond * v
    // see eq (1) in Iwasa Pomiankowski 1994 Evolution 48: 853-867
    double s = t + tprime * v;

    // using open - ended preferences
    return(exp(parms.a * p * s));
} // end Simulation::calculate_attractiveness

// produce a bunch of offspring
void Simulation::offspring_production()
{

}
