/*
 * Basic GL Window modified from the example here
 * http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
 * adapted to use NGL
 */
#include "NGLScene.h"
#include <QKeyEvent>
#include <QApplication>
#include <memory>
#include <iostream>
#include <ngl/NGLInit.h>
#include <ngl/ShaderLib.h>

constexpr float gridSize=1.5;
constexpr int steps=24;
constexpr auto shaderProgram = "Grid";

NGLScene::NGLScene()
{
  setTitle("Simulation");
}

NGLScene::~NGLScene()
{
  glDeleteBuffers(1,&m_vboPointer);
}

void NGLScene::initializeGL()
{
  ngl::NGLInit::instance();
  glewInit();
  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
  initGridShaders();
  makeGrid();
}

void NGLScene::paintGL()
{
  glViewport(0,0,m_win.width,m_win.height);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawGrid();
  update();
}

void NGLScene::initGridShaders()
{
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();

  constexpr auto vertexShader  = "GridVertex";
  constexpr auto fragShader    = "GridFragment";

  shader->createShaderProgram( shaderProgram );

  shader->attachShader( vertexShader, ngl::ShaderType::VERTEX );
  shader->attachShader( fragShader, ngl::ShaderType::FRAGMENT );

  shader->loadShaderSource( vertexShader, "shaders/GridVertex.glsl" );
  shader->loadShaderSource( fragShader, "shaders/GridFragment.glsl" );

  shader->compileShader( vertexShader );
  shader->compileShader( fragShader );

  shader->attachShaderToProgram( shaderProgram, vertexShader );
  shader->attachShaderToProgram( shaderProgram, fragShader );
  
  shader->linkProgramObject( shaderProgram );

  ( *shader )[ shaderProgram ]->use();

  ngl::Vec3 from{ 0.0f, 2.0f, 2.0f };
  ngl::Vec3 to{ 0.0f, 0.0f, 0.0f };
  ngl::Vec3 up{ 0.0f, 1.0f, 0.0f };

  m_view=ngl::lookAt(from,to,up);

  std::cout << "Successfully initialized grid shader." << std::endl;
}

void NGLScene::loadMatricesToShader()
{
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  shader->use("Grid");
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

void  NGLScene::makeGrid()
{
  GLfloat _size=gridSize;
  size_t _steps = steps;
  m_vboSize= (_steps+2)*12;
  std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);

  float step=_size/static_cast<float>(_steps);
	float s2=_size/2.0f;
	float v=-s2;

  for(size_t i=0; i<=_steps; ++i)
	{

    m_gridVBO.push_back({-s2, v, 0.0f});
    m_gridVBO.push_back({s2, v, 0.0f});
    m_gridVBO.push_back({v, s2, 0.0f});
    m_gridVBO.push_back({v, -s2, 0.0f});

		v+=step;
	}

  size_t size = m_gridVBO.size();

  m_gridVAO = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
  m_gridVAO->bind();

  m_gridVAO->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_gridVBO[0].m_x));
  m_gridVAO->setNumIndices(size);
  m_gridVAO->setVertexAttributePointer(0,3,GL_FLOAT,0,0);
  glPointSize(30);
  m_gridVAO->setMode(GL_POINTS);

  m_gridVAO->unbind();
}

void NGLScene::drawGrid()
{
  // std::cout << "Drawing the grid." << std::endl;

  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  
  loadMatricesToShader();

  shader->use(shaderProgram);

  m_gridVAO->bind();
  m_gridVAO->draw();
  m_gridVAO->unbind();
}

void NGLScene::timerEvent(QTimerEvent *)
{
  update();
}

void NGLScene::keyPressEvent(QKeyEvent *_event)
{
  switch (_event->key())
  {
   case Qt::Key_Escape : QApplication::exit(EXIT_SUCCESS); break;
  }
}

void NGLScene::resizeGL(int _w, int _h)
{
  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
}
