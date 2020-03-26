#version 410 core
// this demo is based on code from here https://learnopengl.com/#!PBR/Lighting
/// @brief the vertex passed in
layout (location = 0) in vec3 inVert;
/// @brief the normal passed in
layout (location = 1) in vec3 inNormal;

out vec4 normal;
out vec4 colourNormal;

layout(std140) uniform TransformUBO
{
  mat4 MVP;
  mat4 normalMatrix;
  mat4 M;
}transforms;

void main()
{
    gl_Position = transforms.MVP*vec4(inVert,1.0);
    normal = transforms.MVP*vec4(inNormal, 0.0);
    colourNormal = vec4(inNormal, 1.0f);
}
