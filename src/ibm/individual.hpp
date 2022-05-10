#ifndef _INDIVIDUAL_HPP_
#define _INDIVIDUAL_HPP_

class Individual
{
    public:
        // ornament and preference loci, [2] reflects diploidy
        //
        // ornament locus reflecting local adaptation
        double t_loc[2] = {0.0,0.0};


        // preference locus reflecting choice 
        // based on local adaptation
        double p_loc[2] = {0.0,0.0};


        // ornament locus reflecting male-male competition
        double t_comp[2] = {0.0,0.0};
        
        // preference locus reflecting choice 
        // based on male-male competition
        double p_comp[2] = {0.0,0.0};

        double v_env[2] = {0.0,0.0};

        Individual();
        
        Individual(
                double const p_loc_init, 
                double const t_loc_init, 
                double const p_comp_init, 
                double const t_comp_init,
                double const v_env);

        Individual(Individual const &other);

        void operator=(Individual const &other);


};


#endif
