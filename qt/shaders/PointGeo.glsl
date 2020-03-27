#version 330 core
layout(points) in;
layout(triangle_strip, max_vertices = 10) out;

in vec4 normal[];
in vec4 colourNormal[];
in vec4 perpNormal[];
in vec4 z[];

uniform float normalSize;

out vec4 colour;

void main()
{
    colour=colourNormal[0];

    vec4 inPos = gl_in[0].gl_Position;
    vec4 inNormal = normal[0];
    vec4 extrudeAmount = inNormal*normalSize*0.45;
    vec4 centrePoint = inPos+extrudeAmount;
    float perpSize = normalSize*0.1;
    float size = perpSize;

    // Front
    gl_Position = inPos + (perpNormal[0] + z[0])*perpSize;
    EmitVertex();
    gl_Position = centrePoint;
    EmitVertex();
    gl_Position = inPos + (-perpNormal[0]+z[0])*perpSize;
    EmitVertex();

    // Right
    gl_Position = inPos - (perpNormal[0] + z[0])*perpSize;
    EmitVertex();
    gl_Position = centrePoint;
    EmitVertex();

    // Back
    gl_Position = inPos + (perpNormal[0] - z[0])*perpSize;
    EmitVertex();
    gl_Position = centrePoint;
    EmitVertex();

    // Left
    gl_Position = inPos + (perpNormal[0] + z[0])*perpSize;
    EmitVertex();
    gl_Position = centrePoint;
    EmitVertex();

    EndPrimitive();
}

