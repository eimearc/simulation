#ifndef OPENGLWINDOW_H_
#define OPENGLWINDOW_H_

#include <GL/glew.h>
#include <QOpenGLWindow>


class OpenGLWindow : public QOpenGLWindow
{
  // need to tell Qt to run the MOC
  Q_OBJECT
  public:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief ctor for OpenGL window must set the surface type to OpenGL
    /// @param [in] parent the parent window to the class
    //----------------------------------------------------------------------------------------------------------------------
    explicit OpenGLWindow();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief dtor, remember to remove the device once finished
    //----------------------------------------------------------------------------------------------------------------------

    ~OpenGLWindow();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief render method called every update
    //----------------------------------------------------------------------------------------------------------------------
    void paintGL();
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief pure virtual initialize method we override in our base class to do our drawing
    /// this is only called one time, just after we have a valid GL context use this to init any global GL elements
    //----------------------------------------------------------------------------------------------------------------------
   void initializeGL();
  private:
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief event called by the timer to allow use to re-draw / animate
    //----------------------------------------------------------------------------------------------------------------------
    void timerEvent(QTimerEvent *);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief process key events
    //----------------------------------------------------------------------------------------------------------------------
    void keyPressEvent(QKeyEvent *);
    //----------------------------------------------------------------------------------------------------------------------
    /// @brief resize event
    //----------------------------------------------------------------------------------------------------------------------
    void resizeGL(int _w, int _h);
    //----------------------------------------------------------------------------------------------------------------------
     /// @brief a simple draw grid function
     /// @param[in] _size the size of the grid (width and height)
     /// @param[in] _step sxstep the spacing between grid points
     /// @param[out] o_dataSize the size of the buffer allocated
     /// @returns a pointer to the allocated VBO
     void  makeGrid(GLfloat _size, size_t _steps);
     /// @brief a pointer to our VBO data
     GLuint m_vboPointer=0;
     /// @brief store the size of the vbo data
     GLint m_vboSize=0;
     //----------------------------------------------------------------------------------------------------------------------
     /// @brief window width from resize event
     //----------------------------------------------------------------------------------------------------------------------
     int m_width;
     //----------------------------------------------------------------------------------------------------------------------
     /// @brief window height from resize event
     //----------------------------------------------------------------------------------------------------------------------
     int m_height;

};

#endif
