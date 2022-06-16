#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

// parameters used in the simulation with some default values
struct Parameters
{
    // initial value of the preference locus
    // based on local adaptation
    double p_loc_init = 0.0;
    // initial value of the ornament locus
    // based on local adaptation
    double t_init = 0.0;
    
    double t_conddep_init = 0.0;

    // initial value of the preference locus
    // based on male-male competition
    double p_comp_init = 0.0;
    // initial value of the ornament locus
    // based on male-male competition
    double t_comp_init = 0.0;

    // initial value of the local adaptation locus
    double v_env_init = 0.5;

    // efficacy of mate choice
    // see Iwasa et al 1991 Evolution 45: 1431
    double a = 1.0;

    int clutch_size = 0; //Original = 10

    // strength of survival selection against ornamentation
    // see Iwasa et al 1991 Evolution 45: 1431
    double ct = 0.01;
    double cp = 0.01;

    double mu_p = 0.01;
    double mu_t = 0.01;
    double mu_v = 0.01;
    double mu_t_conddep = 0.01;
    double sdmu = 0.01;

    double vstrength = 0.01;

    // number of males and females
    // per patch
    unsigned nm = 10; //Original =10
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

    // data output interval
    long unsigned output_interval = 10;

    // spatial variance in carrying capacity
    double sigma_k = 0.5;

    // the prefix of the file name
    std::string base_name = "sim_sexsel_space";
};

#endif
