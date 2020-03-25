#include "NGLScene.h"
#include <ngl/NGLInit.h>
#include <ngl/NGLStream.h>
#include <ngl/ShaderLib.h>
#include <ngl/VAOPrimitives.h>
#include <QGuiApplication>
#include <QMouseEvent>


NGLScene::NGLScene()
{
  setTitle( "Simulation");
}


NGLScene::~NGLScene()
{
  std::cout << "Shutting down NGL, removing VAO's and Shaders\n";
}


void NGLScene::resizeGL( int _w, int _h )
{
  m_projection=ngl::perspective( 45.0f, static_cast<float>( _w ) / _h, 0.1f, 200.0f );

  m_win.width  = static_cast<int>( _w * devicePixelRatio() );
  m_win.height = static_cast<int>( _h * devicePixelRatio() );
}
constexpr auto shaderProgram = "PBR";


void NGLScene::initializeGL()
{
  ngl::NGLInit::instance();
  glClearColor( 0.4f, 0.4f, 0.4f, 1.0f );
  glEnable( GL_DEPTH_TEST );
  glEnable( GL_MULTISAMPLE );
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();

  constexpr auto vertexShader  = "PBRVertex";
  constexpr auto fragShader    = "PBRFragment";

  shader->createShaderProgram( shaderProgram );
  shader->attachShader( vertexShader, ngl::ShaderType::VERTEX );
  shader->attachShader( fragShader, ngl::ShaderType::FRAGMENT );
  shader->loadShaderSource( vertexShader, "shaders/PBRVertex.glsl" );
  shader->loadShaderSource( fragShader, "shaders/PBRFragment.glsl" );
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
  shader->setUniform( "camPos", from );
  m_lightPos.set( 0.0, 2.0f, 2.0f ,1.0f);
  shader->setUniform("lightPosition",m_lightPos.toVec3());
  shader->setUniform("lightColor",400.0f,400.0f,400.0f);
  shader->setUniform("exposure",2.2f);
  shader->setUniform("albedo",0.950f, 0.71f, 0.29f);

  shader->setUniform("metallic",1.02f);
  shader->setUniform("roughness",0.38f);
  shader->setUniform("ao",0.2f);
  ngl::VAOPrimitives::instance()->createTrianglePlane("floor",20,20,1,1,ngl::Vec3::up());
  shader->printRegisteredUniforms(shaderProgram);
  shader->use(ngl::nglCheckerShader);
  shader->setUniform("lightDiffuse",1.0f,1.0f,1.0f,1.0f);
  shader->setUniform("checkOn",true);
  shader->setUniform("lightPos",m_lightPos.toVec3());
  shader->setUniform("colour1",0.9f,0.9f,0.9f,1.0f);
  shader->setUniform("colour2",0.6f,0.6f,0.6f,1.0f);
  shader->setUniform("checkSize",60.0f);
  shader->printRegisteredUniforms(ngl::nglCheckerShader);

}


void NGLScene::loadMatricesToShader()
{
  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  shader->use("PBR");
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
   if(m_transformLight)
   {
     shader->setUniform("lightPosition",(m_mouseGlobalTX*m_lightPos).toVec3());

   }
}

void NGLScene::paintGL()
{
  glViewport( 0, 0, m_win.width, m_win.height );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  ngl::ShaderLib* shader = ngl::ShaderLib::instance();
  ( *shader )[ "PBR" ]->use();

  ngl::Mat4 rotX;
  ngl::Mat4 rotY;
  rotX.rotateX( m_win.spinXFace );
  rotY.rotateY( m_win.spinYFace );
  m_mouseGlobalTX = rotX * rotY;
  m_mouseGlobalTX.m_m[ 3 ][ 0 ] = m_modelPos.m_x;
  m_mouseGlobalTX.m_m[ 3 ][ 1 ] = m_modelPos.m_y;
  m_mouseGlobalTX.m_m[ 3 ][ 2 ] = m_modelPos.m_z;

  ngl::VAOPrimitives* prim = ngl::VAOPrimitives::instance();
  loadMatricesToShader();
  prim->draw( "teapot" );
  shader->use(ngl::nglCheckerShader);
  ngl::Mat4 tx;
  tx.translate(0.0f,-0.45f,0.0f);
  ngl::Mat4 MVP=m_projection*m_view*m_mouseGlobalTX*tx;
  ngl::Mat3 normalMatrix=m_view*m_mouseGlobalTX;
  normalMatrix.inverse().transpose();
  shader->setUniform("MVP",MVP);
  shader->setUniform("normalMatrix",normalMatrix);
  if(m_transformLight)
  {
    shader->setUniform("lightPosition",(m_mouseGlobalTX*m_lightPos).toVec3());

  }
  prim->draw("floor");

}


void NGLScene::keyPressEvent( QKeyEvent* _event )
{
  switch ( _event->key() )
  {
    case Qt::Key_Escape:
      QGuiApplication::exit( EXIT_SUCCESS );
      break;
#ifndef USINGIOS_
    case Qt::Key_W:
      glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
      break;
    case Qt::Key_S:
      glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
      break;
#endif
    case Qt::Key_F:
      showFullScreen();
      break;
    case Qt::Key_N:
      showNormal();
      break;
    case Qt::Key_Space :
      m_win.spinXFace=0;
      m_win.spinYFace=0;
      m_modelPos.set(ngl::Vec3::zero());
    break;

    case Qt::Key_L :
    m_transformLight^=true;
    break;
    default:
      break;
  }
  update();
}
