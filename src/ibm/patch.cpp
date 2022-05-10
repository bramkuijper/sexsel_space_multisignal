#include <vector>
#include "patch.hpp"
#include "simulation.hpp"
#include "individual.hpp"

Patch::Patch()
{}

Patch::Patch(int const nm
        ,Individual const &a_male
        ,int const nf
        ,Individual const &a_female
        ,double const envt
        ) :
    males(nm,a_male)
    ,females(nf,a_female)
    ,env{env}
{}
