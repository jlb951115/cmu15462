#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  transformation = svg_2_screen;

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  this->w = this->target_w * sample_rate;
  this->h = this->target_w * sample_rate;
  this->sample_buffer.resize(4 * this->w * this->h, 0);

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;
  this->w = width;
  this->h = height;
  this->sample_rate = 1;
  this->sample_buffer = std::vector<unsigned char>(4 * this->w * this->h);

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack
  Matrix3x3 trans_save = transformation;
  transformation = transformation * element->transform;

  switch(element->type) {
    case POINT:
      draw_point(static_cast<Point&>(*element));
      break;
    case LINE:
      draw_line(static_cast<Line&>(*element));
      break;
    case POLYLINE:
      draw_polyline(static_cast<Polyline&>(*element));
      break;
    case RECT:
      draw_rect(static_cast<Rect&>(*element));
      break;
    case POLYGON:
      draw_polygon(static_cast<Polygon&>(*element));
      break;
    case ELLIPSE:
      draw_ellipse(static_cast<Ellipse&>(*element));
      break;
    case IMAGE:
      draw_image(static_cast<Image&>(*element));
      break;
    case GROUP:
      draw_group(static_cast<Group&>(*element));
      break;
    default:
      break;
  }

  transformation = trans_save;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  float rb = color.r;
  float gb = color.g;
  float bb = color.b;
  float ab = color.a;
  float ra = render_target[4 * (sx + sy * target_w)] / 255.0f;
  float ga = render_target[4 * (sx + sy * target_w) + 1] / 255.0f;
  float ba = render_target[4 * (sx + sy * target_w) + 2] / 255.0f;
  float aa = render_target[4 * (sx + sy * target_w) + 3] / 255.0f;
  float rc = rb * ab + (1 - ab) * aa * ra;
  float gc = gb * ab + (1 - ab) * aa * ga;
  float bc = bb * ab + (1 - ab) * aa * ba;
  float ac = ab + (1 - ab) * aa;
  render_target[4 * (sx + sy * target_w)] = (uint8_t) (rc * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (gc * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (bc * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (ac * 255);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

  // Task 2: 
  // Implement line rasterization
  if (this->sample_rate == 1)
  {
    if (x0 == x1)
    {
      int y_begin = (int) floor(y0);
      int y_end = (int) floor(y1);
      if (y_begin > y_end)
      {
        int temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      for (int y = y_begin; y <= y_end; y++)
        rasterize_point(x0, y, color);
      return;
    }
    if (y0 == y1)
    {
      int x_begin = (int) floor(x0);
      int x_end = (int) floor(x1);
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
      }
      for (int x = x_begin; x <= x_end; x++)
        rasterize_point(x, y0, color);
      return;
    }
    float k = (y1 - y0) / (x1 - x0);
    if (k > 0 && k <= 1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int x = x_begin; x < x_end; x++)
      {
        rasterize_point(x, y_begin, color);
        if (eps + k < 0.5)
          eps = eps + k;
        else
        {
          eps = eps + k - 1;
          y_begin++;
        }
      }
    }
    else if (k > 1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (y_begin > y_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int y = y_begin; y < y_end; y++)
      {
        rasterize_point(x_begin, y, color);
        if (eps + 1.0 / k < 0.5)
          eps = eps + 1.0 / k;
        else
        {
          eps = eps + 1.0 / k - 1;
          x_begin++;
        }
      }
    }
    else if (k < 0 && k >= -1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int x = x_begin; x < x_end; x++)
      {
        rasterize_point(x, y_begin, color);
        if (eps + k > -0.5)
          eps = eps + k;
        else
        {
          eps = eps + k + 1;
          y_begin--;
        }
      }
    }
    else
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (y_begin > y_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int y = y_begin; y < y_end; y++)
      {
        rasterize_point(x_begin, y, color);
        if (eps + 1.0 / k > -0.5)
          eps = eps + 1.0 / k;
        else
        {
          eps = eps + 1.0 / k + 1;
          x_begin--;
        }
      }
    }
  }
  else
  {
     if (x0 == x1)
    {
      int y_begin = (int) floor(y0);
      int y_end = (int) floor(y1);
      if (y_begin > y_end)
      {
        int temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      for (int y = y_begin; y <= y_end; y++)
        fill_pixel(x0, y, color);
      return;
    }
    if (y0 == y1)
    {
      int x_begin = (int) floor(x0);
      int x_end = (int) floor(x1);
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
      }
      for (int x = x_begin; x <= x_end; x++)
        fill_pixel(x, y0, color);
      return;
    }
    float k = (y1 - y0) / (x1 - x0);
    if (k > 0 && k <= 1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int x = x_begin; x < x_end; x++)
      {
        fill_pixel(x, y_begin, color);
        if (eps + k < 0.5)
          eps = eps + k;
        else
        {
          eps = eps + k - 1;
          y_begin++;
        }
      }
    }
    else if (k > 1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (y_begin > y_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int y = y_begin; y < y_end; y++)
      {
        fill_pixel(x_begin, y, color);
        if (eps + 1.0 / k < 0.5)
          eps = eps + 1.0 / k;
        else
        {
          eps = eps + 1.0 / k - 1;
          x_begin++;
        }
      }
    }
    else if (k < 0 && k >= -1)
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (x_begin > x_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int x = x_begin; x < x_end; x++)
      {
        fill_pixel(x, y_begin, color);
        if (eps + k > -0.5)
          eps = eps + k;
        else
        {
          eps = eps + k + 1;
          y_begin--;
        }
      }
    }
    else
    {
      int x_begin = (int) floor(x0);       
      int x_end = (int) floor(x1);       
      int y_begin = (int) floor(y0);       
      int y_end = (int) floor(y1);       
      if (y_begin > y_end)
      {
        int temp = x_begin;
        x_begin = x_end;
        x_end = temp;
        temp = y_begin;
        y_begin = y_end;
        y_end = temp;
      }
      float eps = 0.0;
      for (int y = y_begin; y < y_end; y++)
      {
        fill_pixel(x_begin, y, color);
        if (eps + 1.0 / k > -0.5)
          eps = eps + 1.0 / k;
        else
        {
          eps = eps + 1.0 / k + 1;
          x_begin--;
        }
      }
    }
  }
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 3: 
  // Implement triangle rasterization
  if (this->sample_rate == 1)
  {
  float x_min = std::min(x0, std::min(x1, x2));
  float x_max = std::max(x0, std::max(x1, x2));
  float y_max = std::max(y0, std::max(y1, y2));
  float y_min = std::min(y0, std::min(y1, y2));
  for (int x = (int)floor(x_min); x <= (int)floor(x_max + 1); x++)
  {
    for (int y = (int)floor(y_min); y <= (int) floor(y_max + 1); y++)
    {
      float z1 = (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0);
      float z2 = (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1);
      float z3 = (x0 - x2) * (y - y2) - (x - x2) * (y0 - y2);
      if ((z1 >= 0 && z2 >= 0 && z3 >= 0) || (z1 <= 0 && z2 <= 0 && z3 <= 0))
        rasterize_point(x, y, color);
    }
  }
  } else
  {
  x0 = x0 * sample_rate;
  x1 = x1 * sample_rate;
  x2 = x2 * sample_rate;
  y0 = y0 * sample_rate;
  y1 = y1 * sample_rate;
  y2 = y2 * sample_rate;
  float x_min = std::min(x0, std::min(x1, x2));
  float x_max = std::max(x0, std::max(x1, x2));
  float y_max = std::max(y0, std::max(y1, y2));
  float y_min = std::min(y0, std::min(y1, y2));
  for (int x = (int)floor(x_min); x <= (int)floor(x_max + 1); x++)
  {
    for (int y = (int)floor(y_min); y <= (int) floor(y_max + 1); y++)
    {
      float z1 = (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0);
      float z2 = (x2 - x1) * (y - y1) - (x - x1) * (y2 - y1);
      float z3 = (x0 - x2) * (y - y2) - (x - x2) * (y0 - y2);
      if ((z1 >= 0 && z2 >= 0 && z3 >= 0) || (z1 <= 0 && z2 <= 0 && z3 <= 0))
        fill_sample(x, y, color);
    }
  }
  }

}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization
  float u_scale = tex.width / (x1 - x0);
  float v_scale = tex.height / (y1 - y0);
  for (float x = x0 ; x <= x1; x++)
  {
    for (float y = y0; y <= y1; y++)
    {
      float u = abs(x - x0) / (x1 - x0);
      float v = abs(y - y0) / (y1 - y0);
      rasterize_point(x, y, this->sampler->sample_trilinear(tex, u, v, u_scale, v_scale));
    }
  }
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  if (this->sample_rate == 1)
    return;
  else
  {
    for (size_t i = 0; i < target_w; i++)
    {
      for (size_t j = 0; j < target_h; j++)
      {
        float r = 0.0;
        float g = 0.0;
        float b = 0.0;
        float a = 0.0;
        for (size_t k = 0; k < sample_rate; k++)
        {
          for (size_t m = 0; m < sample_rate; m++)
          {
            size_t x = i * sample_rate + k;
            size_t y = j * sample_rate + m;
            r = r + sample_buffer[4 * (x + y * w)]/ 255.0;
            g = g + sample_buffer[4 * (x + y * w) + 1]/ 255.0;
            b = b + sample_buffer[4 * (x + y * w) + 2]/ 255.0;
            a = a + sample_buffer[4 * (x + y * w) + 3]/ 255.0;
          }
        }
        r = r / sample_rate / sample_rate;
        g = g / sample_rate / sample_rate;
        b = b / sample_rate / sample_rate;
        a = a / sample_rate / sample_rate;
        render_target[4 * (i + j * target_w)] = (uint8_t) (r * 255.0);
        render_target[4 * (i + j * target_w) + 1] = (uint8_t) (g * 255.0);
        render_target[4 * (i + j * target_w) + 2] = (uint8_t) (b * 255.0);
        render_target[4 * (i + j * target_w) + 3] = (uint8_t) (a * 255.0);
      }
    }
  }
  memset(&sample_buffer[0], 0, sample_buffer.size() * sizeof(sample_buffer[0]));

}

void SoftwareRendererImp::fill_sample(int x, int y, const Color& c)
{
  if (x < 0 || x >= this->w) return;
  if (y < 0 || y >= this->h) return;
  sample_buffer[4 * (x + y * w)] = (uint8_t) (c.r * 255.0);
  sample_buffer[4 * (x + y * w) + 1] = (uint8_t) (c.g * 255.0);
  sample_buffer[4 * (x + y * w) + 2] = (uint8_t) (c.b * 255.0);
  sample_buffer[4 * (x + y * w) + 3] = (uint8_t) (c.a * 255.0);
}

void SoftwareRendererImp::fill_pixel(int x, int y, const Color& c)
{
  for (size_t i = 0; i < sample_rate; i++)
    for (size_t j = 0; j < sample_rate; j++)
      fill_sample(x * sample_rate + i, y * sample_rate + j, c);
}

} // namespace CMU462
