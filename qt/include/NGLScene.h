#pragma once

#include <GL/glew.h>
#include <ngl/Text.h>
#include <ngl/Vec3.h>
#include <ngl/Vec4.h>
#include <ngl/Mat4.h>
#include "WindowParams.h"
#include <QOpenGLWindow>
#include <ngl/AbstractVAO.h>
#include <string>
#include <Point.h>
#include <Grid.h>
#include <VectorField.h>

class NGLScene : public QOpenGLWindow
{
  Q_OBJECT
public:
    NGLScene();
    ~NGLScene()=default;
    void initializeGL() override;
    void paintGL() override;
    void resizeGL(int _w, int _h) override;

private:
    Grid m_grid;
    VectorField m_vectorField;
    float m_stepSize;

    WinParams m_win;
    ngl::Mat4 m_mouseGlobalTX;
    ngl::Vec3 m_modelPos;
    ngl::Mat4 m_view;
    ngl::Mat4 m_projection;
    bool m_transformLight=false;
    ngl::Vec4 m_lightPos;
    bool m_drawGrid=false;

    void drawVectorField();
    void initShaders();
    void loadMatricesToShader(const std::string &_shaderName);
    void timerEvent(QTimerEvent *) override;

    // Mouse controls.
    void keyPressEvent(QKeyEvent *_event) override;
    void mouseMoveEvent(QMouseEvent* _event) override;
    void mousePressEvent(QMouseEvent* _event) override;
    void mouseReleaseEvent(QMouseEvent* _event) override;
    void wheelEvent(QWheelEvent* _event) override;
};
