#include "Grid.h"

#include <GLFW/glfw3.h>
#include <gtest/gtest.h>
#include <ngl/NGLStream.h>
#include <iostream>
#include <ngl/NGLInit.h>
#include "Util.h"

using namespace ::testing;

TEST(Grid, startCoords)
{
    setup();

    size_t width = 2;
    size_t height = 3;
    size_t depth = 4;
    float size=1.0f;

    Grid grid(width, height, depth, size);

    ngl::Vec3 expected(0.5f, -0.5f, -0.5f);
    ngl::Vec3 got;

    grid.startCoords(got);

    EXPECT_EQ(expected, got);
}

TEST(Grid, lines)
{
    setup();

    size_t width = 2;
    size_t height = 2;
    size_t depth = 2;
    float size = 1.0f;

    Grid grid(width, height, depth, size);

    for (auto p : grid.m_vbo)
    {
//        std::cout << p << std::endl;
    }
}
