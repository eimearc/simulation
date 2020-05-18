#include <Util.h>

#include <ngl/ShaderLib.h>

void initShader(const std::string &_name, bool _geo)
{
    ngl::ShaderLib *shader = ngl::ShaderLib::instance();

    const std::string vertexShader = _name + "Vertex";
    const std::string geoShader = _name + "Geo";
    const std::string fragShader = _name + "Fragment";

    shader->createShaderProgram(_name);

    shader->attachShader(vertexShader, ngl::ShaderType::VERTEX);
    if (_geo) shader->attachShader(geoShader, ngl::ShaderType::GEOMETRY);
    shader->attachShader(fragShader, ngl::ShaderType::FRAGMENT);

    shader->loadShaderSource(vertexShader, "shaders/"+vertexShader+".glsl");
    if (_geo) shader->loadShaderSource(geoShader, "shaders/"+geoShader+".glsl");
    shader->loadShaderSource(fragShader, "shaders/"+fragShader+".glsl");

    shader->compileShader(vertexShader);
    if (_geo) shader->compileShader(geoShader);
    shader->compileShader(fragShader);

    shader->attachShaderToProgram(_name, vertexShader);
    if (_geo) shader->attachShaderToProgram(_name, geoShader);
    shader->attachShaderToProgram(_name, fragShader);

    shader->linkProgramObject(_name);
}

//void loadMatricesToShader(const std::string &_shaderName, ngl::Mat4 _mouseGlobal, ngl::Mat4 _view, ngl::Mat4 _projection)
//{
//  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
//  shader->use(_shaderName);

//  struct transform
//  {
//    ngl::Mat4 MVP;
//    ngl::Mat4 normalMatrix;
//    ngl::Mat4 M;
//  };

//  transform t;
//  t.M=_view*_mouseGlobal;
//  t.MVP=_projection*t.M;
//  t.normalMatrix=t.M;
//  t.normalMatrix.inverse().transpose();
//  shader->setUniformBuffer("TransformUBO",sizeof(transform),&t.MVP.m_00);
//}


