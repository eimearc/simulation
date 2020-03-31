#version 330 core
layout(points) in;
layout(triangle_strip, max_vertices = 10) out;

in vec4 normal[];
in vec4 colourNormal[];
in vec4 perpNormalU[];
in vec4 perpNormalV[];

uniform float stepSize;

out vec4 colour;

float sq(float x)
{
    return x*x;
}

float length(vec3 v)
{
    return sqrt(sq(v.x) + sq(v.y) + sq(v.y));
}

void main()
{
    colour=colourNormal[0];

    vec4 inPos = gl_in[0].gl_Position;
    float magnitude = length(normal[0]);
    magnitude = clamp(magnitude, stepSize/10, stepSize);
    vec4 inNormal = normalize(normal[0]);
    vec4 normalU = normalize(perpNormalU[0]);
    vec4 normalV = normalize(perpNormalV[0]);

    vec4 extrudeAmount = inNormal*magnitude*0.5;
    vec4 frontCentrePoint = inPos + extrudeAmount;
    vec4 backCentrePoint = inPos - extrudeAmount;

    float perpSize = stepSize*0.05;

    // Front
    gl_Position = backCentrePoint + (normalU + normalV)*perpSize;
    EmitVertex();
    gl_Position = frontCentrePoint;
    EmitVertex();
    gl_Position = backCentrePoint + (-normalU + normalV)*perpSize;
    EmitVertex();

    // Right
    gl_Position = backCentrePoint - (normalU + normalV)*perpSize;
    EmitVertex();
    gl_Position = frontCentrePoint;
    EmitVertex();

    // Back
    gl_Position = backCentrePoint + (normalU - normalV)*perpSize;
    EmitVertex();
    gl_Position = frontCentrePoint;
    EmitVertex();

    // Left
    gl_Position = backCentrePoint + (normalU + normalV)*perpSize;
    EmitVertex();
    gl_Position = frontCentrePoint;
    EmitVertex();

    EndPrimitive();
}

