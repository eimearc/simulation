#include "Point.h"

#include <gtest/gtest.h>
#include <iostream>
#include <Util.h>

using namespace ::testing;

TEST(Point, ctor)
{
    setup();

    ngl::Vec3 position(0,0,0);
    ngl::Vec3 direction(0,1,0);
    ngl::Vec3 velocity(0,1,0);

    Point point(position, direction, velocity);

    EXPECT_EQ(point.position(), position);
    EXPECT_EQ(point.velocity(), velocity);
    EXPECT_EQ(point.direction(), direction);
}

TEST(Point, update)
{
    setup();

    ngl::Vec3 position(0,0,0);
    ngl::Vec3 direction(0,1,0);
    ngl::Vec3 velocity(0,1,0);

    Point point(position, direction, velocity);
    point.update();

    EXPECT_EQ(point.position(), position);
    EXPECT_EQ(point.velocity(), velocity);
    EXPECT_EQ(point.direction(), ngl::Vec3(-0.00999983,0.99995,0));
}
