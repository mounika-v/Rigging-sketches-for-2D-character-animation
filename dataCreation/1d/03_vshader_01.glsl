#version 330

attribute vec4 vPosition;
attribute vec4 vColor;
out vec4 color;
uniform mat4 uModelViewMatrix;

void main (void) 
{
  gl_Position = vPosition;
  gl_PointSize = 10.0;
  color = vColor;
}
