#version 330 core
// This code is based on code from here https://learnopengl.com/#!PBR/Lighting
layout (location =0) out vec4 fragColour;

in vec4 colourOut;

void main()
{
    fragColour = vec4(1.0f, 0.0f, 0.0f, 1.0);
    fragColour = colourOut;
}
