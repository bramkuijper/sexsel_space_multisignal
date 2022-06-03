#include <random>
#include "individual.hpp"

Individual::Individual()
{}

// non-default constructor that builds an
// Individual() object based on a set of
// arguments
Individual::Individual(
        double const p_loc_init, // value of the preference locus based on local adaptation
        double const t_init, // value of the ornament locus (baseline)
        double const t_conddep_init, // value of the ornament locus based on local adaptation
        double const p_comp_init, // value of the preference locus based on competition
        double const t_comp_init, // value of the ornament locus based on competition
        double const v_env // value of the local adaptation locus
        ) :
    p_loc{p_loc_init,p_loc_init}
    ,t{t_init,t_init}
    ,t_conddep{t_conddep_init,t_conddep_init}
    ,p_comp{p_comp_init,p_comp_init}
    ,t_comp{t_comp_init,t_comp_init}
    ,v_env{v_env,v_env}
{}

// copy constructor
Individual::Individual(
        Individual const &other) :
    p_loc{other.p_loc[0],other.p_loc[1]}
    ,t{other.t[0],other.t[1]}
    ,t_conddep{other.t_conddep[0],other.t_conddep[1]}
    ,p_comp{other.p_comp[0],other.p_comp[1]}
    ,t_comp{other.t_comp[0],other.t_comp[1]}
    ,v_env{other.v_env[0],other.v_env[1]}
{
}

// offspring production constructor
Individual::Individual(
    Individual &mom
    ,Individual dad
    ,std::mt19937 &rng
    ,std::bernoulli_distribution &allele_sampler) :
    p_loc{mom.p_loc[allele_sampler(rng)],
        dad.p_loc[allele_sampler(rng)]}

    ,t{mom.t[allele_sampler(rng)],
        dad.t[allele_sampler(rng)]}

    ,t_conddep{mom.t_conddep[allele_sampler(rng)],
        dad.t_conddep[allele_sampler(rng)]}

    ,p_comp{mom.p_comp[allele_sampler(rng)],
        dad.p_comp[allele_sampler(rng)]}

    ,t_comp{mom.t_comp[allele_sampler(rng)],
        dad.t_comp[allele_sampler(rng)]}

    ,v_env{mom.v_env[allele_sampler(rng)],
        dad.v_env[allele_sampler(rng)]}
{}

//assignment operator
void Individual::operator=(Individual const &other)
{
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        p_loc[allele_idx] = other.p_loc[allele_idx];
        t[allele_idx] = other.t[allele_idx];
        t_conddep[allele_idx] = other.t_conddep[allele_idx]; 
        p_comp[allele_idx] = other.p_comp[allele_idx];
        t_comp[allele_idx] = other.t_comp[allele_idx];
        v_env[allele_idx] = other.v_env[allele_idx];
    }
}


// mutate loci
void Individual::mutate(
        std::mt19937 &rng_r // random number generator
        ,Parameters &params) // parameter object
{
    std::uniform_real_distribution<double> uniform_dist{0.0,1.0};
    std::normal_distribution<> normal{0.0,params.sdmu};

    // go through each allele and mutate it with a certain probability
    // if mutation applies, then change allelic value by adding a number
    // drawn from a normal distribution
    for (int allele_idx = 0; allele_idx < 2; ++allele_idx)
    {
        if (uniform_dist(rng_r) < params.mu_p)
        {
            p_loc[allele_idx] += normal(rng_r);
        }

        if (uniform_dist(rng_r) < params.mu_t)
        {
            t[allele_idx] += normal(rng_r);
        }
        
        if (uniform_dist(rng_r) < params.mu_t_conddep)
        {
            t_conddep[allele_idx] += normal(rng_r);
        }
        
        if (uniform_dist(rng_r) < params.mu_v)
        {
            v_env[allele_idx] += normal(rng_r);
        }
    } // end for allele_idx;

} // end Individual::mutate()
