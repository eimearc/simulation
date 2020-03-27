#version 330 core
layout(points) in;
layout(line_strip, max_vertices = 2) out;

in vec4 normal[];
in vec4 colourNormal[];

uniform float normalSize;

out vec4 colour;

void main()
{
    colour=colourNormal[0];

    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[0].gl_Position + (normal[0]*normalSize);
    EmitVertex();

    EndPrimitive();
}

