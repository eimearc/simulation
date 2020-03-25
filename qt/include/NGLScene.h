#ifndef NGLSCENE_H
#define NGLSCENE_H

#include <GL/glew.h>
#include <ngl/Text.h>
#include <ngl/Vec3.h>
#include <ngl/Vec4.h>
#include <ngl/Mat4.h>
#include "WindowParams.h"
#include <QOpenGLWindow>
#include <ngl/AbstractVAO.h>

class NGLScene : public QOpenGLWindow
{
  Q_OBJECT
public:
  NGLScene();
  ~NGLScene() override;
  void initializeGL() override;
  void paintGL() override;
  void resizeGL(int _w, int _h) override;

private:
  GLuint m_vboPointer=0;
  GLint m_vboSize=0;

  std::unique_ptr<ngl::AbstractVAO> m_gridVAO;
  std::vector<ngl::Vec3> m_gridVBO;

  WinParams m_win;
  ngl::Mat4 m_mouseGlobalTX;
  ngl::Vec3 m_modelPos;
  ngl::Mat4 m_view;
  ngl::Mat4 m_projection;
  bool m_transformLight=false;
  ngl::Vec4 m_lightPos;

  void makeGrid();
  void drawGrid();
  void initGridShaders();
  void loadMatricesToShader();
  void keyPressEvent(QKeyEvent *_event) override;
  void timerEvent(QTimerEvent *);
};

#endif
