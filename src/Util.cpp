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


