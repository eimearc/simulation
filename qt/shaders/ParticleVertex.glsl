#version 410 core
// this demo is based on code from here https://learnopengl.com/#!PBR/Lighting
/// @brief the vertex passed in
layout (location = 0) in vec2 inVert;

const float numParticles = 1000.0f;

out vec4 color;

layout(std140) uniform TransformUBO
{
  mat4 MVP;
  mat4 normalMatrix;
  mat4 M;
}transforms;

void main()
{
    gl_Position = transforms.MVP*vec4(inVert,0.0f,1.0);

    float r = ((gl_VertexID & 0x0000000F) >> 0) / 16.0f;
    float g = ((gl_VertexID & 0x000000F0) >> 4) / 16.0f;
    float b = ((gl_VertexID & 0x00000F00) >> 8) / 16.0f;

    color = vec4(r,g,b,1);
}
