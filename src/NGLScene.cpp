#include "NGLScene.h"

#include "Util.h"
#include <QKeyEvent>
#include <QApplication>
#include <memory>
#include <iostream>
#include <ngl/NGLInit.h>
#include <ngl/ShaderLib.h>
#include <ngl/NGLStream.h>
#include <cstdlib>
#include <gflags/gflags.h>

DEFINE_int32(resolution, 50, "Resolution of the MAC grid.");

constexpr float GRID_SIZE=1;
constexpr auto GRID_SHADER = "Grid";
constexpr auto POINT_SHADER = "Point";
constexpr auto PARTICLE_SHADER = "Particle";

NGLScene::NGLScene()
{
  setTitle("Simulation");
}

void NGLScene::initializeGL()
{
  ngl::NGLInit::instance();
  glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
  glEnable(GL_DEPTH_TEST);

  initShaders();

  m_stepSize = GRID_SIZE/static_cast<float>(FLAGS_resolution);
  m_grid = Grid(FLAGS_resolution, FLAGS_resolution, 0, GRID_SIZE);
  m_vectorField = VectorField(FLAGS_resolution, FLAGS_resolution, 0, GRID_SIZE);
  m_macGrid = MAC(FLAGS_resolution);
}

void NGLScene::paintGL()
{
  glViewport(0,0,m_win.width,m_win.height);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  ngl::Mat4 rotX;
  ngl::Mat4 rotY;
  rotX.rotateX( m_win.spinXFace );
  rotY.rotateY( m_win.spinYFace );
  m_mouseGlobalTX = rotX * rotY;
  m_mouseGlobalTX.m_m[ 3 ][ 0 ] = m_modelPos.m_x;
  m_mouseGlobalTX.m_m[ 3 ][ 1 ] = m_modelPos.m_y;
  m_mouseGlobalTX.m_m[ 3 ][ 2 ] = m_modelPos.m_z;

  drawMACGrid();

  if (m_drawGrid)
  {
    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    loadMatricesToShader(GRID_SHADER);
    shader->setUniform("colour", ngl::Vec3(0.1,0.4,0.1));
    m_grid.drawInner();
  }
}

void NGLScene::drawMACGrid()
{
    loadMatricesToShader(PARTICLE_SHADER);
    loadMatricesToShader(GRID_SHADER);
    m_macGrid.update();
    m_macGrid.draw();
    update();
}

void NGLScene::drawVectorField()
{
    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    shader->use(POINT_SHADER);
    shader->setUniform("stepSize", m_stepSize);
    loadMatricesToShader(POINT_SHADER);

    m_vectorField.update();
    m_vectorField.draw();
    update();
}

void NGLScene::initShaders()
{
  initShader(GRID_SHADER);
  initShader(POINT_SHADER, true);
  initShader(PARTICLE_SHADER);

  ngl::Vec3 from{ 0.0f, 0.0f, 2.0f };
  ngl::Vec3 to{ 0.0f, 0.0f, 0.0f };
  ngl::Vec3 up{ 0.0f, 1.0f, 0.0f };
  m_view=ngl::lookAt(from,to,up);
}

void NGLScene::loadMatricesToShader(const std::string &_shaderName)
{
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  shader->use(_shaderName);

  struct transform
  {
    ngl::Mat4 MVP;
    ngl::Mat4 normalMatrix;
    ngl::Mat4 M;
  };

  transform t;
  t.M=m_view*m_mouseGlobalTX;
  t.MVP=m_projection*t.M;
  t.normalMatrix=t.M;
  t.normalMatrix.inverse().transpose();
  shader->setUniformBuffer("TransformUBO",sizeof(transform),&t.MVP.m_00);
}

void NGLScene::timerEvent(QTimerEvent *)
{
  update();
}

void NGLScene::resizeGL(int _w, int _h)
{
  m_projection=ngl::perspective( 45.0f, static_cast<float>( _w ) / _h, 0.1f, 200.0f );
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
}
