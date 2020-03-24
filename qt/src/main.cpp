/****************************************************************************
basic OpenGL demo modified from http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
****************************************************************************/
#include "OpenGLWindow.h"

#include <QtGui/QGuiApplication>
#include <iostream>



int main(int argc, char **argv)
{
  QGuiApplication app(argc, argv);
  // create an OpenGL format specifier
  QSurfaceFormat format;
  // set the number of samples for multisampling
  // will need to enable glEnable(GL_MULTISAMPLE); once we have a context
  format.setSamples(4);
  format.setMajorVersion(4);
  format.setMinorVersion(1);
  // now we are going to set to Compat Profile OpenGL so we can use and old Immediate mode GL
  format.setProfile(QSurfaceFormat::CompatibilityProfile);
  // now set the depth buffer to 24 bits
  format.setDepthBufferSize(24);
  QSurfaceFormat::setDefaultFormat(format);
  // now we are going to create our scene window
  OpenGLWindow window;
  // we can now query the version to see if it worked
  std::cout<<"Profile is "<<format.majorVersion()<<" "<<format.minorVersion()<<"\n";
  // set the window size
  window.resize(1024, 720);
  // and finally show
  window.show();

  return app.exec();
}



