#version 330 core
layout(points) in;
layout(triangle_strip, max_vertices = 3) out;

in vec4 normal[];
in vec4 colourNormal[];
in vec4 perpNormal[];

uniform float normalSize;

out vec4 colour;

void main()
{
    colour=colourNormal[0];

    vec4 inPos = gl_in[0].gl_Position;
    vec4 inNormal = normal[0];
    vec4 extrudeAmount = inNormal*normalSize*0.4;

    gl_Position = inPos + extrudeAmount;
    EmitVertex();
    gl_Position = inPos + perpNormal[0]*0.1;
    EmitVertex();
    gl_Position = inPos - perpNormal[0]*0.1;
    EmitVertex();

    EndPrimitive();
}

