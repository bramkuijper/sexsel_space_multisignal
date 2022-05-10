#ifndef _PATCH_HPP_
#define _PATCH_HPP_

#include <vector>
#include "individual.hpp"

class Patch
{
    public:

        // the population of males and females
        std::vector<Individual> males;
        std::vector<Individual> females;

        double env;

        // default constructor, which sets everything to a default of 0
        Patch();

        // patch-constructor where one can modulate numbers of males and females
        Patch(
                int const nm // number of males in this patch
                ,Individual const &a_male
                ,int const nf // number of females in this patch
                ,Individual const &a_female
                ,double const env // environmental variable
                );
};

#endif
