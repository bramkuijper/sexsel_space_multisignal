#include <gtest/gtest.h>
#include "individual.hpp"

TEST(Mate_Choice_Space_Ind_Test, TestAssignment)
{
    Individual ind1;

    Individual ind2;

    ind2.p_loc[0] = 2.0;
    ind2.t_loc[0] = 3.0;
    ind2.t_comp[0] = 4.0;
    ind2.p_comp[0] = 5.0;

    ind1 = ind2;

    EXPECT_EQ(ind1.p_loc[0], ind1.p_loc[0]);
    EXPECT_EQ(ind1.t_loc[0], ind1.t_loc[0]);
    EXPECT_EQ(ind1.p_comp[0], ind1.p_comp[0]);
    EXPECT_EQ(ind1.p_comp[0], ind1.p_comp[0]);
}
