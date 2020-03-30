#version 330 core
// This code is based on code from here https://learnopengl.com/#!PBR/Lighting
layout (location =0) out vec4 fragColour;

in vec4 colour;

void main()
{
    fragColour = colour;
}
