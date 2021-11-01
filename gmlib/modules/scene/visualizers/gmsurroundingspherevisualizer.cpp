/**********************************************************************************
**
** Copyright (C) 1994 Narvik University College
** Contact: GMlib Online Portal at http://episteme.hin.no
**
** This file is part of the Geometric Modeling Library, GMlib.
**
** GMlib is free software: you can redistribute it and/or modify
** it under the terms of the GNU Lesser General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** GMlib is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU Lesser General Public License for more details.
**
** You should have received a copy of the GNU Lesser General Public License
** along with GMlib.  If not, see <http://www.gnu.org/licenses/>.
**
**********************************************************************************/



// local
#include "gmsurroundingspherevisualizer.h"

//GMlib
#include <scene/camera/gmcamera.h>
#include <scene/render/gmdefaultrenderer.h>



namespace GMlib {

  SurroundingSphereVisualizer::SurroundingSphereVisualizer(int m1, int m2): _m1(m1), _m2(m2), _no_indices(2*(2*_m1-1)*_m2), _line_width(1.0f)  {

    _prog.acquire("color");
    _vbo.create();
    _ibo.create();

    _color = GMcolor::white();
  }




  void SurroundingSphereVisualizer::setSphere(const Sphere<float,3>& s) {

    _sphere = s;

    // Uppdate VBO and grid lines
    _fillVBO();
    _makeLines();
  }



  void SurroundingSphereVisualizer::render(const SceneObject *obj, const DefaultRenderer* renderer) const{

      GL_CHECK(::glLineWidth(_line_width));

      _prog.bind(); {
        _prog.uniform( "u_mvpmat", obj->getModelViewProjectionMatrix(renderer->getCamera()));
        _prog.uniform( "u_color", _color) ;
        _prog.uniform( "u_selected", false );

        GL::AttributeLocation vert_loc = _prog.getAttributeLocation( "in_vertex" );

        _vbo.bind();
          _vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<const GLvoid*>(0x0) );

          _ibo.bind();
            glDrawElements( GL_LINES, _no_indices, GL_UNSIGNED_SHORT, reinterpret_cast<const GLvoid*>(0x0) );
          _ibo.unbind();

          _vbo.disable( vert_loc );
        _vbo.unbind();
      } _prog.unbind();
  }




  void SurroundingSphereVisualizer::update() {

      GL::GLVertex *ptr = _vbo.mapBuffer<GL::GLVertex>();
      DVector<GL::GLVertex> dp(ptr, 2+(_m1-1)*_m2, 0);
      _fill(dp);
      _vbo.unmapBuffer();
  }




  void SurroundingSphereVisualizer::_fillVBO() {

    // Fill the vertice buffer
    DVector<GL::GLVertex> dp(2+(_m1-1)*_m2);
    _fill(dp);
    const GLsizeiptr data_size =  sizeof(GL::GLVertex) * static_cast<unsigned int>(dp.getDim());
    _vbo.bufferData( data_size, dp.getPtr(), GL_DYNAMIC_DRAW );
  }




  void SurroundingSphereVisualizer::_fill(DVector<GL::GLVertex>& dp) {

    float x = _sphere.getPos()[0];
    float y = _sphere.getPos()[1];
    float z = _sphere.getPos()[2];
    float r = _sphere.getValue();

    // Compute top of the sphere
    dp[0] = {x, y, z+r};
    // Compute bottom of the sphere
    dp[1] = {x, y, z-r};

    // Compute steps in both parameter directions.
    const double du = M_PI/_m1;
    const double dv = M_2PI/_m2;

    // Compute body points on the sphere
    for( int i=1, k=2; i < _m1; i++ ) {
        const double u  = M_PI_2 - i*du;
        const double cu = cos(u);
        const float  su = r*float(sin(u));
        for( int j = 0; j < _m2; j++ ) {
            const double v = j*dv;
            const float a = float(cu*cos(v));
            const float b = float(cu*sin(v));
            dp[k++] = {x+r*a , y+r*b, z+su};
        }
    }
  }




  void SurroundingSphereVisualizer::_makeLines()     // Fill IBO
  {
      // Create the indice buffer
      std::vector<GLushort> indices(_no_indices);

      unsigned int k=0;
      for(int i=0; i < _m2; i++) {
          indices[k++] = GLushort(0);
          indices[k++] = GLushort(i+2);
          for(int j=0; j< _m1-2; j++) {
              indices[k++] = GLushort(2+i+j*_m2);
              indices[k++] = GLushort(2+i+(j+1)*_m2);
          }
          indices[k++] = GLushort(2+i+(_m1-2)*_m2);
          indices[k++] = GLushort(1);
      }
      for(int j=0; j< _m1-1; j++) {
          for(int i=0; i < _m2-1; i++) {
              indices[k++] = GLushort(2+i+j*_m2);
              indices[k++] = GLushort(3+i+j*_m2);
          }
          indices[k++] = GLushort(1+(j+1)*_m2);
          indices[k++] = GLushort(2+j*_m2);
      }
      _ibo.bufferData( GLushort(_no_indices) * sizeof(GLushort), indices.data(), GL_DYNAMIC_DRAW );
  }

}// END namespace GMlib
