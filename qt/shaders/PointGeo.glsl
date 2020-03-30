#version 330 core
layout(points) in;
layout(triangle_strip, max_vertices = 10) out;

in vec4 normal[];
in vec4 colourNormal[];
in vec4 perpNormalU[];
in vec4 perpNormalV[];

uniform float normalSize;

out vec4 colour;

void main()
{
    colour=colourNormal[0];

    vec4 inPos = gl_in[0].gl_Position;
    vec4 inNormal = normal[0];

    float magnitude = sqrt(inNormal.x*inNormal.x + inNormal.y*inNormal.y + inNormal.z*inNormal.z);
    inNormal = normalize(inNormal);

    vec4 extrudeAmount = inNormal*magnitude*normalSize*0.2;
    vec4 frontCentrePoint = inPos+extrudeAmount;
    vec4 backCentrePoint = inPos - extrudeAmount;
    float perpSize = normalSize;
    perpSize = clamp(perpSize, 0.01, 0.02);

    vec4 normalU = normalize(perpNormalU[0]);
    vec4 normalV = normalize(perpNormalV[0]);

    // Front
    gl_Position = backCentrePoint + (normalU + normalV)*perpSize;
    EmitVertex();
    gl_Position = frontCentrePoint;
    EmitVertex();
    gl_Position = backCentrePoint + (-normalU + normalV)*perpSize;
    EmitVertex();

    // Right
    gl_Position = backCentrePoint- (normalU + normalV)*perpSize;
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

