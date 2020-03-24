#include <iostream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

GLFWwindow *window;

void initWindow()
{
}

void mainLoop()
{
    while(!glfwWindowShouldClose(window))
    {
    }
}

void cleanup()
{
    glfwDestroyWindow(window);
    glfwTerminate();
}

void run()
{
    initWindow();
    mainLoop();
    cleanup();
}

int main()
{
    std::cout << "Hello, world!" << std::endl;

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL 
    window = glfwCreateWindow(500, 500, "Simulation", nullptr, nullptr);
    if (window == NULL)
    {
        throw std::runtime_error("Failed to open GLFW window.");
    }
    glfwMakeContextCurrent(window);
    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

    do{

        // Draw nothing, see you in tutorial 2 !

        // Swap buffers
        glfwSwapBuffers(window);
        glfwPollEvents();

    } // Check if the ESC key was pressed or the window was closed
    while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
        glfwWindowShouldClose(window) == 0 );
}