#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        // ornament and preference loci, [2] reflects diploidy
        //
        // ornament locus reflecting local adaptation
        double t[2] = {0.0,0.0};


        // preference locus reflecting choice 
        // based on local adaptation
        double p_loc[2] = {0.0,0.0};

        // signaling how strongly ornament size should scale with
        // local adaptation. This trait is analagous to t' in
        // eq (1) in Iwasa Pomiankowski 1994 Evolution 48: 853-867
        // whereas t above is analagous to t in that eq (1)
        double t_conddep[2] = {0.0,0.0};


        // ornament locus reflecting male-male competition
        // (not used)
        double t_comp[2] = {0.0,0.0};
        
        // preference locus reflecting choice 
        // based on male-male competition
        double p_comp[2] = {0.0,0.0};

        double v_env[2] = {0.0,0.0};

        Individual();

        // initialization constructor
        // if you want nonzero values to initialize stuff
        Individual(
                double const p_loc_init, 
                double const t_init, 
                double const t_conddep_init, 
                double const p_comp_init, 
                double const t_comp_init,
                double const v_env);

        // copy constructor
        Individual(Individual const &other);

        // offspring production constructor
        Individual(Individual &mother
                ,Individual father
                ,std::mt19937 &rng
                ,std::bernoulli_distribution &allele_sampler);


        void operator=(Individual const &other);

        // mutation function
        void mutate(std::mt19937 &random_number_generator
                ,Parameters &params);

};


#endif
