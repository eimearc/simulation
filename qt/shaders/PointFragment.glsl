#version 330 core
// This code is based on code from here https://learnopengl.com/#!PBR/Lighting
layout (location =0) out vec4 fragColour;

in vec4 colour;

void main()
{
//    fragColour = vec4(0.0f, 0.5f, 0.6f, 1.0);
    fragColour = colour;
}
