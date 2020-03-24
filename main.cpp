#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

GLFWwindow *window;

void initWindow()
{
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
}

void mainLoop()
{
    float vertices[] = {
        -0.5f, -0.5f, 0.0f,
        0.5f, -0.5f, 0.0f,
        0.0f,  0.5f, 0.0f
    };

    unsigned int VBO;
    // glGenBuffers(1, &VBO);

    while(!glfwWindowShouldClose(window))
    {
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glfwSwapBuffers(window);
        glfwPollEvents();
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

    run();
}