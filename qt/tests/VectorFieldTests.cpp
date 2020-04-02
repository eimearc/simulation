#include "VectorField.h"

#include <gtest/gtest.h>
#include <iostream>
#include "Util.h"

using namespace ::testing;

TEST(VectorField, ctor)
{
    setup();

    size_t width = 2;
    size_t height = 2;
    size_t depth = 2;
    float size = 1.0f;

    VectorField vectorField(width, height, depth, size);
    Grid grid(width, height, depth, size);

    ngl::Vec3 expected(-0.25f, -0.25, -0.25f);
    ngl::Vec3 got;

    vectorField.startCoords(got);

    EXPECT_EQ(expected, got);
    EXPECT_EQ(grid, vectorField.m_grid);
}
