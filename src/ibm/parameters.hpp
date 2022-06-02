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
    double t_loc_init = 0.0;
    // initial value of the preference locus
    // based on male-male competition
    double p_comp_init = 0.0;
    // initial value of the ornament locus
    // based on male-male competition
    double t_comp_init = 0.0;

    // initial value of the local adaptation locus
    double v_env_init = 0.5;

    // strength of survival selection against ornamentation
    double ct = 0.01;
    double cp = 0.01;

    double vstrength = 0.01;

    // number of males and females
    // per patch
    unsigned nm = 10;
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
