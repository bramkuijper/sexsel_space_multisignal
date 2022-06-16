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

        // have offspring disperse across space and have them
        // replace fathers and mothers
        dispersal_and_replacement();

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
        << "nf;" << parms.nf << std::endl
        << "a;" << parms.a << std::endl
        << "ct;" << parms.ct << std::endl
        << "cp;" << parms.cp << std::endl
        << "clutch_size;" << parms.clutch_size << std::endl

        << "mu_p;" << parms.mu_p << std::endl
        << "mu_t;" << parms.mu_t << std::endl
        << "mu_t_conddep;" << parms.mu_t_conddep << std::endl
        << "mu_v;" << parms.mu_v << std::endl
        << "sdmu;" << parms.sdmu << std::endl

        << "npatches;" << parms.Npatches << std::endl
        << "d;" << parms.d << std::endl
        << "k0;" << parms.k0 << std::endl
        << "a;" << parms.a << std::endl
        << "sigma_k;" << parms.sigma_k << std::endl;
} // write_parameters()

// write the stats to a file
void Simulation::write_data()
{
    // set up variables to take stats
    double mean_t = 0;
    double var_t = 0;

    double mean_t_conddep = 0;
    double var_t_conddep = 0;

    double mean_p = 0;
    double var_p = 0;

    double mean_v = 0;
    double var_v = 0;
    
    double mean_s = 0;
    double var_s = 0;

    double p,t,tpr,v,s;
    int n = 0;

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

            mean_p += p;
            var_p += p * p;

            mean_t += t;
            var_t += t * t;
            
            mean_t_conddep += tpr;
            var_t_conddep += tpr * tpr;

            mean_v += v;
            var_v += v * v;

            mean_s += s;
            var_s += s * s;
        } // end iterate over females
        
        for (std::vector<Individual>::iterator male_iter = patch_iter->males.begin(); 
                male_iter != patch_iter->males.end(); 
                ++male_iter)
        {
            p = male_iter->p_loc[0] + male_iter->p_loc[1];
            t = male_iter->t[0] + male_iter->t[1];
            tpr = male_iter->t_conddep[0] + male_iter->t_conddep[1];
            v = male_iter->v_env[0] + male_iter->v_env[1];
            s = t + tpr * v;

            mean_p += p;
            var_p += p * p;

            mean_t += t;
            var_t += t * t;
            
            mean_t_conddep += tpr;
            var_t_conddep += tpr * tpr;

            mean_v += v;
            var_v += v * v;

            mean_s += s;
            var_s += s * s;
        } // end iterate over males
    
        n += patch_iter->males.size() + patch_iter->females.size();
    } // end iterate over patches


    mean_p /= n;
    mean_t /= n;
    mean_t_conddep /= n;
    mean_s /= n;
    mean_v /= n;

    // variance is 
    var_p = var_p / n - mean_p * mean_p;
    var_t = var_t / n - mean_t * mean_t;
    var_t_conddep = var_t_conddep / n - mean_t_conddep * mean_t_conddep;
    var_s = var_s / n - mean_s * mean_s;
    var_v  = var_v / n - mean_v * mean_v;

    output_file << time_step << ";"
        << mean_survivors_f << ";"
        << mean_survivors_m << ";"


        << mean_p << ";"
        << mean_t << ";"
        << mean_t_conddep << ";"
        << mean_s << ";"
        << mean_v << ";"

        << var_p << ";"
        << var_t << ";"
        << var_t_conddep << ";"
        << var_s << ";"
        << var_v << ";";

    output_file << std::endl;
} // end Simulation::write_data()


// initialize the data files used to output the simulation results
void Simulation::initialize_output_file()
{
    output_file << "generation;mean_survivors_f;mean_survivors_m;mean_p;mean_t;mean_t_conddep;mean_s;mean_v;var_p;var_t;var_t_conddep;var_s;var_v;" << std::endl;
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
                    ,parms.t_init
                    ,parms.t_conddep_init
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
    double t = male_i.t[0] + male_i.t[1];

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
// we are not going to implement this yet, first focus on 
// female choice and carrying capacity only
void Simulation::male_male_competition()
{}

// females choose males and make kids
void Simulation::female_choice()
{
    // make a distribution to sample either allele 1 vs allele 2 
    // when inheriting
    std::bernoulli_distribution diploid_sampler{0.5};

    int male_chosen;

    // reset juvenile counter
    njuveniles = 0;

    // go through all the patches and perform choice
//    for (std::vector<Patch>::iterator patch_iter = metapop.begin(); 
//            patch_iter != metapop.end(); 
//            ++patch_iter)
    for (int patch_idx = 0;
            patch_idx < metapop.size(); 
            ++patch_idx)
    {
        // remove any previous juveniles
        // as we are going to add new ones
        metapop[patch_idx].juveniles.clear();

        // patch extinct, hence this patch will not produce offspring
        // continue to the next patch
        if (metapop[patch_idx].females.size() == 0 || 
                metapop[patch_idx].males.size() == 0)
        {
            continue;
        }

        for (std::vector<Individual>::iterator female_iter = 
                metapop[patch_idx].females.begin(); 
                female_iter != metapop[patch_idx].females.end(); 
                ++female_iter)
        {
            std::vector<double> attractiveness_values;

            // loop through males and make a distribution 
            // of all the attractiveness values
            for (std::vector<Individual>::iterator male_iter = 
                    metapop[patch_idx].males.begin(); 
                    male_iter != metapop[patch_idx].males.end(); 
                    ++male_iter)
            {
                // calculate attractiveness and add this to a vector
                attractiveness_values.push_back(
                        calculate_attractiveness(*female_iter, *male_iter)
                        );
            } // end iterate over all males

            assert(attractiveness_values.size() == metapop[patch_idx].males.size());
            // set up a probability distribution of all the attractiveness
            // values and choose from them
            // males with larger attractiveness values 
            // are more likely to be chosen
            std::discrete_distribution<int> choose_male_dist(
                    attractiveness_values.begin()
                    ,attractiveness_values.end());

            // now sample a male from this distribution
            male_chosen = choose_male_dist(rng_r);

            // make sure this gives a sensible value
            assert(male_chosen >= 0);
            assert(male_chosen < metapop[patch_idx].males.size());

            for (int egg_idx = 0; egg_idx < parms.clutch_size; ++egg_idx)
            {
                // make offspring from mom and dad
                // see the individual.hpp's offspring constructor
                Individual kid(*female_iter
                        ,metapop[patch_idx].males[male_chosen]
                        ,rng_r
                        ,diploid_sampler
                        );
                
                kid.mutate(rng_r,parms);

                // now have this pair produce offspring
                metapop[patch_idx].juveniles.push_back(kid);

                ++njuveniles;
            }
        } // end iterate over all females
    } // end for patch
} // end Simulation::female_choice()


// calculate attractiveness values
double Simulation::calculate_attractiveness(Individual &female
        ,Individual &male)
{
    // express ornament
    double p = female.p_loc[0] + female.p_loc[1];
    double t = male.t[0] + male.t[1];
    double tprime = male.t_conddep[0] + male.t_conddep[1];
    double v = male.v_env[0] + male.v_env[1];

    // express male ornament as a condition-dependent trait
    // hence we assume s = t + tcond * v
    // see eq (1) in Iwasa Pomiankowski 1994 Evolution 48: 853-867
    double s = t + tprime * v;

    // using open - ended preferences
    return(exp(parms.a * p * s));
} // end Simulation::calculate_attractiveness


void Simulation::dispersal_and_replacement()
{
    // go through all the patches and make a distribution of 
    // juveniles that we will sample when a dispersal event happens
    //
    // patches with lots of juveniles
    std::vector<int> dispersal_vector;

    for (std::vector<Patch>::iterator patch_iter = metapop.begin(); 
            patch_iter != metapop.end(); 
            ++patch_iter)
    {
        dispersal_vector.push_back(patch_iter->juveniles.size());
    }

    // make a dispersal probabilitty distribution
    // so that patches with more juveniles are more likely to be 
    // sampled when sampling immigrants
    std::discrete_distribution<int> patch_of_origin_sampler{dispersal_vector.begin()
        ,dispersal_vector.end()};

    int patch_of_origin_idx, sampled_juvenile_idx;


    // no juveniles to replace adults
    if (njuveniles < 1)
    {
        std::cout << "population extinct in generation " << time_step << ": no juveniles available." << std::endl;
        write_parameters();
        exit(1);
    }

    // ok, now sample juveniles either from local or remote populations
    for (int patch_idx = 0;
            patch_idx < metapop.size(); 
            ++patch_idx)
    {
        // remove current males and females
        metapop[patch_idx].males.clear();
        metapop[patch_idx].females.clear();

        // first sample males from the juvenile population
        for (int male_idx = 0; male_idx < parms.nm; ++male_idx)
        {
            // no locals available?
            // or fate determines we need a dispersing offspring?
            //
            // get immigrant offspring
            if (metapop[patch_idx].juveniles.size() < 1 ||
                    uniform(rng_r) < parms.d)
            {
                do 
                {
                    // sample a remote patch of origin
                    patch_of_origin_idx = patch_of_origin_sampler(rng_r);
                }
                while (metapop[patch_of_origin_idx].juveniles.size() < 1);
                // sample again if there is nothing in this patch
            }
            else
            {
                // patch of origin is local patch
                patch_of_origin_idx = patch_idx;
            }

            // some bounds checking
            assert(patch_of_origin_idx >= 0);
            assert(patch_of_origin_idx <  metapop.size());
            assert(metapop[patch_of_origin_idx].juveniles.size() > 0);

            std::uniform_int_distribution juvenile_sampler(0, 
                (int)metapop[patch_of_origin_idx].juveniles.size() - 1);

            sampled_juvenile_idx = juvenile_sampler(rng_r);
                
            // add juvenile to the stack of males
            metapop[patch_of_origin_idx].males.push_back(
                    metapop[patch_of_origin_idx].juveniles[sampled_juvenile_idx]);
        }
        
        for (int female_idx = 0; female_idx < parms.nf; ++female_idx)
        {
            // no locals available?
            // or fate determines we need a dispersing offspring?
            //
            // get immigrant offspring
            if (metapop[patch_idx].juveniles.size() < 1 ||
                    uniform(rng_r) < parms.d)
            {
                // sample a patch of origin
                patch_of_origin_idx = patch_of_origin_sampler(rng_r);

                // some bounds checking
                assert(patch_of_origin_idx >= 0);
                assert(patch_of_origin_idx <  metapop.size());


                // make a random distribution to sample a juvenile
                // from this remote patch
                std::uniform_int_distribution juvenile_sampler(0, 
                    (int)metapop[patch_of_origin_idx].juveniles.size() - 1);

                // draw a number from this distribution
                sampled_juvenile_idx = juvenile_sampler(rng_r);
            }
            else
            {
                patch_of_origin_idx = patch_idx;

                std::uniform_int_distribution juvenile_sampler(0, 
                    (int)metapop[patch_of_origin_idx].juveniles.size() - 1);

                sampled_juvenile_idx = juvenile_sampler(rng_r);
            }
            
            // add juvenile to the stack of females 
            metapop[patch_of_origin_idx].females.push_back(
                    metapop[patch_of_origin_idx].juveniles[sampled_juvenile_idx]);
        }
        
    }
} // Simulation::dispersal_and_replacement()
