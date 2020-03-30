#include "Point.h"

#include <gtest/gtest.h>
#include <iostream>

using namespace ::testing;

TEST(Point, basic)
{
    std::cout << "Working" << std::endl;
    EXPECT_EQ(1, 1);
}
