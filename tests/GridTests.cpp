#include "Grid.h"

#include <GLFW/glfw3.h>
#include <gtest/gtest.h>
#include <ngl/NGLStream.h>
#include <iostream>
#include <ngl/NGLInit.h>
#include "Util.h"

using namespace ::testing;

TEST(Grid, ctor)
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
    EXPECT_EQ(grid.width(), width);
    EXPECT_EQ(grid.height(), height);
    EXPECT_EQ(grid.depth(), depth);
    EXPECT_EQ(grid.gridSize(), size);
    EXPECT_EQ(grid.stepSize(), size/width);
}

TEST(Grid, lines)
{
    setup();

    size_t width = 2;
    size_t height = 2;
    size_t depth = 2;
    float size = 1.0f;

    Grid grid(width, height, depth, size);

//    const std::vector<ngl::Vec3> &vbo = grid.m_vbo;

    const std::vector<ngl::Vec3> expected = {
        {0.5,-0.5,-0.5},
        {0.5,0.5,-0.5},
        {0,-0.5,-0.5},
        {0,0.5,-0.5},
        {-0.5,-0.5,-0.5},
        {-0.5,0.5,-0.5},
        {0.5,-0.5,0},
        {0.5,0.5,0},
        {0,-0.5,0},
        {0,0.5,0},
        {-0.5,-0.5,0},
        {-0.5,0.5,0},
        {0.5,-0.5,0.5},
        {0.5,0.5,0.5},
        {0,-0.5,0.5},
        {0,0.5,0.5},
        {-0.5,-0.5,0.5},
        {-0.5,0.5,0.5},
        {0.5,-0.5,-0.5},
        {-0.5,-0.5,-0.5},
        {0.5,0,-0.5},
        {-0.5,0,-0.5},
        {0.5,0.5,-0.5},
        {-0.5,0.5,-0.5},
        {0.5,-0.5,0},
        {-0.5,-0.5,0},
        {0.5,0,0},
        {-0.5,0,0},
        {0.5,0.5,0},
        {-0.5,0.5,0},
        {0.5,-0.5,0.5},
        {-0.5,-0.5,0.5},
        {0.5,0,0.5},
        {-0.5,0,0.5},
        {0.5,0.5,0.5},
        {-0.5,0.5,0.5},
        {0.5,-0.5,-0.5},
        {0.5,-0.5,0.5},
        {0.5,0,-0.5},
        {0.5,0,0.5},
        {0.5,0.5,-0.5},
        {0.5,0.5,0.5},
        {0,-0.5,-0.5},
        {0,-0.5,0.5},
        {0,0,-0.5},
        {0,0,0.5},
        {0,0.5,-0.5},
        {0,0.5,0.5},
        {-0.5,-0.5,-0.5},
        {-0.5,-0.5,0.5},
        {-0.5,0,-0.5},
        {-0.5,0,0.5},
        {-0.5,0.5,-0.5},
        {-0.5,0.5,0.5}
    };

//    EXPECT_EQ(vbo.size(), 54);
//    EXPECT_EQ(vbo, expected);
}
