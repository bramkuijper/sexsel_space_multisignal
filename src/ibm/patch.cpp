#include <vector>
#include "patch.hpp"
#include "simulation.hpp"
#include "individual.hpp"

// default patch constructors (i.e., set everything to zero)
Patch::Patch()
{}

//! patch constructor with parameters
//! @param nm Integer denoting the number of males in the local patch
//! @param a_male Individual object reflecting the typical individual that
//!                 makes a male
//! @param nf Integer denoting the number of females in the local patch
//! @param a_female Individual object reflecting the typical individual that
//!                 makes a female
//! @param envt Floating point value reflecting the environmental gradient
//!                 (which is value between 0 and 1)
Patch::Patch(int const nm
        ,Individual const &a_male
        ,int const nf
        ,Individual const &a_female
        ,double const envt
        ) :
    males(nm,a_male)
    ,females(nf,a_female)
    ,coordinate{envt}
{}
