#include "MAC.h"

#include <gtest/gtest.h>

TEST(MAC, ctor)
{
    EXPECT_EQ(1,1);
}

TEST(MACGrid, ctor)
{
    EXPECT_EQ(1,1);

    MAC::Grid grid(3);
    EXPECT_EQ(grid.m_resolution, 3);
}
