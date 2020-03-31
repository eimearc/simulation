#include "Grid.h"

#include <GLFW/glfw3.h>

#include <gtest/gtest.h>
#include <ngl/NGLStream.h>
#include <iostream>
#include <ngl/NGLInit.h>

using namespace ::testing;

void setup()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(100, 100, "Simulation Test", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    ngl::NGLInit::instance();
}

TEST(Grid, startCoords)
{
    setup();

    size_t width = 2;
    size_t height = 3;
    size_t depth = 4;
    float size=1.0f;

    Grid grid(width, height, depth, size);

    ngl::Vec3 expected;
    ngl::Vec3 got;

    grid.startCoords(got);

    EXPECT_EQ(expected, got);
}
