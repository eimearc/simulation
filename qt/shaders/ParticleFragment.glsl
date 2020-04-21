#version 330 core
// This code is based on code from here https://learnopengl.com/#!PBR/Lighting
in vec4 color;

layout (location =0) out vec4 fragColour;

void main()
{
    fragColour = color;
}
