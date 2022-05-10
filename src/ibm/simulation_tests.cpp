#include <gtest/gtest.h>
#include "simulation.hpp"

TEST(Mate_Choice_Space_Sim_Test, TestInitialization)
{
    Parameters pars;
    Simulation the_sim(pars);

    EXPECT_EQ(the_sim.metapop.size(), pars.Npatches);
}
