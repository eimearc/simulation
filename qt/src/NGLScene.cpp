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
constexpr float gridSize=1.5;
constexpr int steps=24;

NGLScene::NGLScene()
{
  setTitle("Qt5 compat profile OpenGL 3.2");
}

NGLScene::~NGLScene()
{
  // now we have finished clear the device
  std::cout<<"deleting buffer\n";
  glDeleteBuffers(1,&m_vboPointer);
}



void NGLScene::initializeGL()
{
  glewInit();

  glClearColor(0.4f, 0.4f, 0.4f, 1.0f);
  makeGrid();
}



void  NGLScene::makeGrid()
{
  GLfloat _size=gridSize;
  size_t _steps = steps;
  m_vboSize= (_steps+2)*12;
  std::unique_ptr<GLfloat []>vertexData( new GLfloat[m_vboSize]);
  int k=-1;
  float step=_size/static_cast<float>(_steps);
	float s2=_size/2.0f;
	float v=-s2;
  for(size_t i=0; i<=_steps; ++i)
	{
		// vertex 1 x,y,z
		vertexData[++k]=-s2; // x
		vertexData[++k]=v; // y
		vertexData[++k]=0.0; // z

		// vertex 2 x,y,z
		vertexData[++k]=s2; // x
		vertexData[++k]=v; // y
		vertexData[++k]=0.0; // z

		// vertex 3 x,y,z
		vertexData[++k]=v;
		vertexData[++k]=s2;
		vertexData[++k]=0.0;

		// vertex 4 x,y,z
		vertexData[++k]=v;
		vertexData[++k]=-s2;
		vertexData[++k]=0.0;
		// now change our step value
		v+=step;
	}

  glGenBuffers(1, &m_vboPointer);
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  glBufferData(GL_ARRAY_BUFFER, m_vboSize*sizeof(GL_FLOAT) ,
               vertexData.get(), GL_STATIC_DRAW);
}

void NGLScene::drawGrid()
{
  glEnableClientState(GL_VERTEX_ARRAY);
  glBindBuffer(GL_ARRAY_BUFFER, m_vboPointer);
  glPointSize(10.0f);
  glVertexPointer(3,GL_FLOAT,0,0);
  glDrawArrays(GL_LINES, 0, m_vboSize/2);
  glDrawArrays(GL_POINTS, 0, m_vboSize/2);
  glDisableClientState(GL_VERTEX_ARRAY);
}

void NGLScene::paintGL()
{
  glViewport(0,0,m_win.width,m_win.height);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  drawGrid();
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
