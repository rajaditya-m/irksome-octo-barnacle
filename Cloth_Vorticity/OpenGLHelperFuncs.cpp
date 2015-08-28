#include <QGLWidget>
#include <QtOpenGL>
#include <gl/glu.h>

inline float getPixelDepth(int pixel_position_x, int pixel_position_y) {
  GLint view_port[4];
  glGetIntegerv(GL_VIEWPORT, view_port);
  float depth;
  glReadPixels(pixel_position_x, view_port[3] - pixel_position_y, 1, 1,
               GL_DEPTH_COMPONENT, GL_FLOAT, &depth);
  return depth;
}