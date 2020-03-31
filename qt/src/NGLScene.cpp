/*
 * Basic GL Window modified from the example here
 * http://qt-project.org/doc/qt-5.0/qtgui/openglwindow.html
 * adapted to use NGL
 */
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

constexpr float gridSize=1.5;
constexpr int steps=25;

constexpr auto gridShader = "Grid";
constexpr auto pointShader = "Point";

NGLScene::NGLScene()
{
  setTitle("Simulation");
  m_2d=false;
}

NGLScene::~NGLScene()
{
}

void NGLScene::initializeGL()
{
  ngl::NGLInit::instance();
  glewInit();
  glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
  glEnable(GL_DEPTH_TEST);

  initShaders();

  m_grid = Grid(gridSize, steps);
  m_vectorField = VectorField(steps, steps, steps, gridSize, steps);

  makePoints();
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

  if (m_drawGrid)
  {
    loadMatricesToShader(gridShader);
    m_grid.draw();
  }

  drawVectorField();
}

void NGLScene::drawVectorField()
{
    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    shader->use(pointShader);
    shader->setUniform("stepSize", m_stepSize);
    loadMatricesToShader(pointShader);
    m_vectorField.update();
    update();
    m_vectorField.draw();
}

void NGLScene::initShaders()
{
  initShader(gridShader);
  initShader(pointShader, true);

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

void NGLScene::getPointStartCoords(ngl::Vec3 &_coords, float &_step)
{
//  getGridStartCoords(_coords, _step);

  _step = gridSize/static_cast<float>(steps);
  m_stepSize = _step;

  m_grid.startCoords(_coords);

  _coords.m_x *= -1;
  _coords += _step/2.0f;
}

void NGLScene::makePoints()
{
    m_pointsVAO = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
    m_points.clear();

    ngl::Vec3 coords;
    float step;
    getPointStartCoords(coords, step);

    ngl::Vec3 position;
    ngl::Vec3 direction;
    ngl::Vec3 velocity = {0.0f,0.1f,0.0f};

    float u = coords.m_x;
    float v = coords.m_y;
    float w = coords.m_z;

    if (m_2d)
    {
        w = 0.0f;
        for (int i = 0; i < steps; ++i)
        {
            for (int j = 0; j < steps; ++j)
            {
                position = ngl::Vec3(u+step*i,v+step*j,w);
                direction = ngl::Vec3(0.0f, 1.0f, 0.0f);
                if (i%3 == 0)
                {
                    direction = ngl::Vec3(0.0f, 0.0f, 1.0f);
                }
                else if(i%3 == 1)
                {
                    direction = ngl::Vec3(1.0f, 0.0f, 0.0f);
                }
                direction.normalize();

                m_points.push_back({position, direction, velocity});
            }
        }
    }

    else
    {
        for (int i = 0; i < steps; ++i)
        {
            for (int j = 0; j < steps; ++j)
            {
                for(int k = 0; k < steps; ++k)
                {
                    position = ngl::Vec3(u+step*i,v+step*j,w+step*k);
                    direction = ngl::Vec3(0.0f, 1.0f, 0.0f);
                    if (k%3 == 0)
                    {
                        direction = ngl::Vec3(0.0f, 0.0f, 1.0f);
                        direction *= 2.0f;
                    }
                    else if(k%3 == 1)
                    {
                        direction = ngl::Vec3(1.0f, 0.0f, 0.0f);
                        direction *= 0.5f;
                    }

                    m_points.push_back({position, direction, velocity});
                }
            }
        }
    }

    updatePointsVBO();
}

void NGLScene::updatePoints()
{
    for (auto &point : m_points)
    {
        point.update();
    }

    updatePointsVBO();
    update();
}

void NGLScene::updatePointsVBO()
{
    m_pointsVBO.clear();

    for (const auto& point : m_points)
    {
        m_pointsVBO.push_back(point.position());
        m_pointsVBO.push_back(point.direction());
    }

    size_t size = m_pointsVBO.size();

    m_pointsVAO->bind();
        m_pointsVAO->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_pointsVBO[0].m_x));
        m_pointsVAO->setNumIndices(size);
        m_pointsVAO->setVertexAttributePointer(0,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),0); // Position.
        m_pointsVAO->setVertexAttributePointer(1,3,GL_FLOAT,2*(GLsizei)sizeof(ngl::Vec3),3); // Direction.
        m_pointsVAO->setMode(GL_POINTS);
    m_pointsVAO->unbind();
}

void NGLScene::drawPoints()
{
    for (auto point : m_points)
    {
        point.update();
    }

    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    shader->use(pointShader);
    shader->setUniform("normalSize", m_stepSize*0.9f);
    loadMatricesToShader(pointShader);

    m_pointsVAO->bind();
    m_pointsVAO->draw();
    m_pointsVAO->unbind();
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
