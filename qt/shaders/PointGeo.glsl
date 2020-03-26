#version 330 core
layout(points) in;
layout(line_strip, max_vertices = 2) out;

in vec4 normal[];

uniform float size;

out vec3 colour;

void main()
{
    colour=vec3(0,1,0);

    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[0].gl_Position + (normal[0]*0.025f);
    EmitVertex();

    EndPrimitive();
}

