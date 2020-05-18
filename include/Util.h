#pragma once

#include <string>
#include <ngl/Mat4.h>

void initShader(const std::string &_name, bool _geo=false);
void loadMatricesToShader(const std::string &_shaderName, ngl::Mat4 _mouseGlobal, ngl::Mat4 _view, ngl::Mat4 _projection);
