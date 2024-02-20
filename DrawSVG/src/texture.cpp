#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);
    for (int w = 0; w < width; w++)
    {
      for (int h = 0; h < height; h++)
      {
        float r = 0.0;
        float g = 0.0;
        float b = 0.0;
        float a = 0.0;
        for (int k = 0; k < 2; k++)
        {
          for (int m = 0; m < 2; m++)
          {
            r += tex.mipmap[startLevel+i-1].texels[4 * (2 * w + k + 2 * w * (2 * h + m))] / 255.0f;
            g += tex.mipmap[startLevel+i-1].texels[4 * (2 * w + k + 2 * w * (2 * h + m)) + 1] / 255.0f;
            b += tex.mipmap[startLevel+i-1].texels[4 * (2 * w + k + 2 * w * (2 * h + m)) + 1] / 255.0f;
            a += tex.mipmap[startLevel+i-1].texels[4 * (2 * w + k + 2 * w * (2 * h + m)) + 1] / 255.0f;
          }
        }
        r = r / 4.0;
        g = g / 4.0;
        b = b / 4.0;
        a = a / 4.0;
        level.texels[4 * (w + h * w)] = (uint8_t) (r * 255);
        level.texels[4 * (w + h * w) + 1] = (uint8_t) (g * 255);
        level.texels[4 * (w + h * w) + 2] = (uint8_t) (b * 255);
        level.texels[4 * (w + h * w) + 3] = (uint8_t) (a * 255);
      }
    }

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

  // Task 6: Implement nearest neighbour interpolation
  
  // return magenta for invalid level
  int x = (int) floor(tex.mipmap[level].width * u);
  int y = (int) floor(tex.mipmap[level].height * v);
  float r = tex.mipmap[level].texels[4 * (x + tex.mipmap[level].width * y)] / 255.0f;
  float g = tex.mipmap[level].texels[4 * (x + tex.mipmap[level].width * y) + 1] / 255.0f;
  float b = tex.mipmap[level].texels[4 * (x + tex.mipmap[level].width * y) + 2] / 255.0f;
  float a = tex.mipmap[level].texels[4 * (x + tex.mipmap[level].width * y) + 3] / 255.0f;
  return Color(r, g, b, a);

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task 6: Implement bilinear filtering

  // return magenta for invalid level
  float xo = tex.mipmap[level].width * u;
  float yo = tex.mipmap[level].height * v;
  int x0 = (int) floor(tex.mipmap[level].width * u);
  int y0 = (int) floor(tex.mipmap[level].height * v);
  int x1 = (int) floor(tex.mipmap[level].width * u + 1);
  int y1 = (int) floor(tex.mipmap[level].height * v + 1);
  if (x0 >= tex.mipmap[level].width - 1 || y0 >= tex.mipmap[level].height - 1)
    return sample_nearest(tex, u, v, level);
  float t = xo - x0;
  float t1 = yo - y0;
  float r = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y0)] / 255.0f;
  float g = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y0) + 1] / 255.0f;
  float b = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y0) + 2] / 255.0f;
  float a = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y0) + 3] / 255.0f;
  float r1 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y0)] / 255.0f;
  float g1 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y0) + 1] / 255.0f;
  float b1 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y0) + 2] / 255.0f;
  float a1 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y0) + 3] / 255.0f;
  float r2 = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y1)] / 255.0f;
  float g2 = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y1) + 1] / 255.0f;
  float b2 = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y1) + 2] / 255.0f;
  float a2 = tex.mipmap[level].texels[4 * (x0 + tex.mipmap[level].width * y1) + 3] / 255.0f;
  float r3 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y1)] / 255.0f;
  float g3 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y1) + 1] / 255.0f;
  float b3 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y1) + 2] / 255.0f;
  float a3 = tex.mipmap[level].texels[4 * (x1 + tex.mipmap[level].width * y1) + 3] / 255.0f;
  float r_inter_y0 = (1 - t) * r  + t * r1;
  float g_inter_y0 = (1 - t) * g  + t * g1;
  float b_inter_y0 = (1 - t) * b  + t * b1;
  float a_inter_y0 = (1 - t) * a  + t * a1;
  float r_inter_y1 = (1 - t) * r2  + t * r3;
  float g_inter_y1 = (1 - t) * g2  + t * g3;
  float b_inter_y1 = (1 - t) * b2  + t * b3;
  float a_inter_y1 = (1 - t) * a2  + t * a3;
  float r_final = (1- t1) * r_inter_y0 + t1 * r_inter_y1;
  float g_final = (1- t1) * g_inter_y0 + t1 * g_inter_y1;
  float b_final = (1- t1) * b_inter_y0 + t1 * b_inter_y1;
  float a_final = (1- t1) * a_inter_y0 + t1 * a_inter_y1;
  return Color(r, g, b, a);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 7: Implement trilinear filtering

  // return magenta for invalid level
  float d = max(log2f(max(u_scale, v_scale)), 0.0f);
  int ld = (int) floor(d);
  if (ld >= tex.mipmap.size())
    return Color(1, 0, 1, 1);
  int ld_high = ld + 1;;
  if (ld_high >= tex.mipmap.size())
    return sample_bilinear(tex, u, v, ld);
  float t1 = d - ld;
  Color c1 = sample_bilinear(tex, u, v, ld);
  Color c2 = sample_bilinear(tex, u, v, ld_high);
  return (1 - t1) * c1 + t1 * c2;
}

} // namespace CMU462
