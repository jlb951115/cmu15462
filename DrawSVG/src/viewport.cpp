#include "viewport.h"

#include "CMU462.h"

namespace CMU462 {

void ViewportImp::set_viewbox( float centerX, float centerY, float vspan ) {

  // Task 5 (part 2): 
  // Set svg coordinate to normalized device coordinate transformation. Your input
  // arguments are defined as normalized SVG canvas coordinates.
  this->centerX = centerX;
  this->centerY = centerY;
  this->vspan = vspan; 
  float top_left_y = this->centerY - this->vspan;
  float top_left_x = this->centerX - this->vspan;
  Matrix3x3 translate = Matrix3x3::identity();
  translate(0,0) = 1.0; translate(0, 1) = 0.0; translate(0, 2) = -top_left_x;
  translate(1,0) = 0.0; translate(1, 1) = 1.0; translate(1, 2) = -top_left_y;
  translate(2,0) = 0.0; translate(2, 1) = 0.0; translate(2, 2) = 1.0;
  Matrix3x3 scale = Matrix3x3::identity();
  scale(0, 0) = 1.0 / (2 * this->vspan); scale(0, 1) = 0.0; scale(0, 2) = 0.0;
  scale(1, 0) = 0.0; scale(1, 1) = 1.0 / (2 * this->vspan); scale(1, 2) = 0.0;
  scale(2, 0) = 0.0; scale(2, 1) = 0; scale(2, 2) = 1.0;
  set_svg_2_norm(scale * translate);
  }

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->centerX -= dx;
  this->centerY -= dy;
  this->vspan *= scale;
  set_viewbox( centerX, centerY, vspan );
}

} // namespace CMU462
