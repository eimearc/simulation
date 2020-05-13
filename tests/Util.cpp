#include "Util.h"

#include <GLFW/glfw3.h>
#include <iostream>
#include <ngl/NGLInit.h>

void setup()
{
    glfwInit();

    GLFWwindow* window = glfwCreateWindow(100, 100, "Simulation Test", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
    }
    glfwMakeContextCurrent(window);
    ngl::NGLInit::instance();
}
