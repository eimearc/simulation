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
#include <ngl/NGLStream.h>

constexpr float gridSize=1.5;
constexpr int steps=25;

constexpr auto gridShader = "Grid";
constexpr auto pointShader = "Point";

NGLScene::NGLScene()
{
  setTitle("Simulation");
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
  makeGrid();
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
    drawGrid();
  }
//  drawTeapot();
  drawPoints();
}

void NGLScene::drawTeapot()
{
  ngl::VAOPrimitives* prim = ngl::VAOPrimitives::instance();
  loadMatricesToShader("PBR");
  prim->draw( "teapot" );
}

void NGLScene::initShader(const std::string &_name, bool _geo=false)
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

void NGLScene::initShaders()
{

  initShader(gridShader);
  initShader(pointShader, true);

  ngl::Vec3 from{ 0.0f, 2.0f, 2.0f };
  ngl::Vec3 to{ 0.0f, 0.0f, 0.0f };
  ngl::Vec3 up{ 0.0f, 1.0f, 0.0f };
  m_view=ngl::lookAt(from,to,up);

  std::cout << "Successfully initialized shaders." << std::endl;
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

void NGLScene::getGridStartCoords(ngl::Vec3 &_coords, float &_step)
{
  _step = gridSize/static_cast<float>(steps);

  _coords.m_x = gridSize/2.0f;
  _coords.m_y = -(_coords.m_x);
  _coords.m_z = -(_coords.m_x);
}

void NGLScene::getPointStartCoords(ngl::Vec3 &_coords, float &_step)
{
  getGridStartCoords(_coords, _step);

  _coords.m_x *= -1;
  _coords += _step/2.0f;
}

void NGLScene::makeGridVBO()
{
    float step;
    ngl::Vec3 pos;
    getGridStartCoords(pos, step);

    float u = pos.m_x;
    float v = pos.m_y;

    for(size_t i=0; i<=steps; ++i)
    {
        float d = -u;
        for (size_t j=0; j<=steps; ++j)
        {
            m_gridVBO.push_back({-u, v, d}); // Left vert
            m_gridVBO.push_back({u, v, d});  // Right vert
            m_gridVBO.push_back({v, u, d});  // Top vert
            m_gridVBO.push_back({v, -u, d}); // Bottom vert

            d+=step;
        }
        v+=step;
    }

    v = -u;
    for(size_t i=0; i<=steps; ++i)
    {
        float d = -u;
        for (size_t j=0; j<=steps; ++j)
        {
            m_gridVBO.push_back({-u, d, v}); // Left vert
            m_gridVBO.push_back({u, d, v});  // Right vert
            m_gridVBO.push_back({v, d, u});  // Top vert
            m_gridVBO.push_back({v, d, -u}); // Bottom vert

            d+=step;
        }
        v+=step;
    }
}

void NGLScene::makeGrid()
{
  makeGridVBO();

  size_t size = m_gridVBO.size();
  m_gridVAO = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
  m_gridVAO->bind();

  m_gridVAO->setData(ngl::SimpleVAO::VertexData(size*sizeof(ngl::Vec3), m_gridVBO[0].m_x));
  m_gridVAO->setNumIndices(size);
  m_gridVAO->setVertexAttributePointer(0,3,GL_FLOAT,0,0);

  glPointSize(5);
  glEnable(GL_LINE_SMOOTH);
  glLineWidth(50);

  m_gridVAO->setMode(GL_LINES);

  m_gridVAO->unbind();
}

void NGLScene::makePoints()
{
    m_pointsVBO.clear();

    ngl::Vec3 coords;
    float step;
    getPointStartCoords(coords, step);

    float u = coords.m_x;
    float v = coords.m_y;
    float w = coords.m_z;

    size_t numPoints = 0;

    for (int i = 0; i < steps; ++i)
    {
        for (int j = 0; j < steps; ++j)
        {
            for(int k = 0; k < steps; ++k)
            {
                numPoints++;
                m_pointsVBO.push_back({u+step*i,v+step*j,w+step*k});
                ngl::Vec3 dir = {0.0f, 1.0f, 1.0f};
                if (k%3 == 0)
                {
                    dir = ngl::Vec3(1.0f, 0.0f, 0.0f);
                }
                dir.normalize();
                m_pointsVBO.push_back(dir);
            }
        }
    }

    std::cout << "Num points: " << numPoints << " Size: " << m_pointsVBO.size() << std::endl;

    size_t size = m_pointsVBO.size();

    m_pointsVAO = ngl::VAOFactory::createVAO("simpleVAO", GL_POINTS);
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
    ngl::ShaderLib* shader = ngl::ShaderLib::instance();
    shader->use(pointShader);
    loadMatricesToShader(pointShader);

    m_pointsVAO->bind();
    m_pointsVAO->draw();
    m_pointsVAO->unbind();
}

void NGLScene::drawGrid()
{
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  shader->use(gridShader);
  loadMatricesToShader(gridShader);

  m_gridVAO->bind();
  m_gridVAO->draw();
  m_gridVAO->unbind();
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
