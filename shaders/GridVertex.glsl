#version 410 core
// this demo is based on code from here https://learnopengl.com/#!PBR/Lighting
/// @brief the vertex passed in
layout (location = 0) in vec3 inVert;
/// @brief the normal passed in
layout (location = 1) in vec3 inNormal;

uniform vec3 colour;

out vec4 colourOut;

layout(std140) uniform TransformUBO
{
  mat4 MVP;
  mat4 normalMatrix;
  mat4 M;
}transforms;

void main()
{
    gl_Position = transforms.MVP*vec4(inVert.x, inVert.y, 0.0 ,1.0);
    colourOut = vec4(colour,1.0f);
}
