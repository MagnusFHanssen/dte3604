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


#include "gmselectorvisualizer.h"

#include "../gmsceneobject.h"
#include "../camera/gmcamera.h"
#include "../render/gmdefaultrenderer.h"

// gmlib
#include "../../opengl/gmopengl.h"


namespace GMlib {


  SelectorVisualizer *SelectorVisualizer::_s_instance = nullptr;


  SelectorVisualizer::SelectorVisualizer( float r, Material mat )
    : _top_bot_verts(0), _mid_strips(0), _mid_strips_verts(0),
      _mat(mat) {

    _prog.acquire("blinn_phong");
    _color_prog.acquire("color");

    _vbo.create();
    _ibo.create();

    makeGeometry( double(r), 8, 8 );
  }

  SelectorVisualizer::SelectorVisualizer(int m1, int m2, float r, Material mat)
    : _top_bot_verts(0), _mid_strips(0), _mid_strips_verts(0),
      _mat(mat) {

    _prog.acquire("blinn_phong");

    _vbo.create();
    _ibo.create();

    makeGeometry( double(r), m1, m2 );
  }

  void SelectorVisualizer::render(const SceneObject* obj, const DefaultRenderer *renderer) const {

    const Camera* cam = renderer->getCamera();
    const HqMatrix<float,3> &mvmat = obj->getModelViewMatrix(cam);
    const HqMatrix<float,3> &pmat  = obj->getProjectionMatrix(cam);

    _prog.bind(); {

      // Model view and projection matrices
      _prog.uniform( "u_mvmat", mvmat );
      _prog.uniform( "u_mvpmat", pmat * mvmat );

      // Lights
      _prog.bindBufferBase( "DirectionalLights",  renderer->getDirectionalLightUBO(), 0 );
      _prog.bindBufferBase( "PointLights",        renderer->getPointLightUBO(), 1 );
      _prog.bindBufferBase( "SpotLights",         renderer->getSpotLightUBO(), 2 );

      // Material data
      _prog.uniform( "u_mat_amb", _mat.getAmb() );
      _prog.uniform( "u_mat_dif", _mat.getDif() );
      _prog.uniform( "u_mat_spc", _mat.getSpc() );
      _prog.uniform( "u_mat_shi", _mat.getShininess() );

      // Shader attribute locations
      GL::AttributeLocation vert_loc = _prog.getAttributeLocation( "in_vertex" );
      GL::AttributeLocation normal_loc = _prog.getAttributeLocation( "in_normal" );

      // Bind and draw
      _vbo.bind();
        _vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid *>(0x0) );
        _vbo.enable( normal_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid *>(sizeof(GL::GLNormal)) );

        // Draw top and bottom caps
        for( int i = 0; i < 2; i++ )
          glDrawArrays( GL_TRIANGLE_FAN, i * _top_bot_verts, _top_bot_verts );

        // Draw body strips
        for( int i = 0; i < _mid_strips; i++ )
          glDrawArrays( GL_TRIANGLE_STRIP, _top_bot_verts*2 + i*_mid_strips_verts, _mid_strips_verts );

        _vbo.disable( normal_loc );
        _vbo.disable( vert_loc );
      _vbo.unbind();

    } _prog.unbind();
  }



  void SelectorVisualizer::renderGeometry( const SceneObject* obj, const Renderer* renderer, const Color& color) const {

    _color_prog.bind(); {
      _color_prog.uniform( "u_mvpmat", obj->getModelViewProjectionMatrix(renderer->getCamera()) );
      _color_prog.uniform( "u_color", color );

      GL::AttributeLocation vertice_loc = _color_prog.getAttributeLocation( "in_vertex" );

      _vbo.bind();
        _vbo.enable( vertice_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid *>(0x0) );

        // Draw top and bottom caps
        for( int i = 0; i < 2; i++ )
          glDrawArrays( GL_TRIANGLE_FAN, i * _top_bot_verts, _top_bot_verts );

        // Draw body strips
        for( int i = 0; i < _mid_strips; i++ )
          glDrawArrays( GL_TRIANGLE_STRIP, _top_bot_verts*2 + i*_mid_strips_verts, _mid_strips_verts );

        _vbo.disable( vertice_loc );
      _vbo.unbind();

    } _color_prog.unbind();
  }



  SelectorVisualizer *SelectorVisualizer::getInstance() {

    if( !_s_instance )
      _s_instance = new SelectorVisualizer;

    return _s_instance;
  }



  void SelectorVisualizer::makeGeometry(double r, int m1, int m2) {

    _top_bot_verts    = m2+2;
    _mid_strips       = m1-2;
    _mid_strips_verts = (m2+1)*2;

    const int no_verts = _top_bot_verts*2 + _mid_strips*_mid_strips_verts;


    DVector<GL::GLVertexNormal> verts(no_verts);
    GL::GLVertexNormal *verts_ptr = verts.getPtr();

    // Compute stride in the spheres u and v parametric direction.
    const double du = M_PI/m1;
    const double dv = M_2PI/m2;


    // Compute top and bottom triangle fans for the sphere
    {
      double u;
      double su, cu, ru;

      const Arrow<float,3> top = Arrow<float,3>(Point<float,3>(0,0,float(r)), Vector<float,3>(0,0,1));
      verts_ptr->x  = top.getPos()(0);
      verts_ptr->y  = top.getPos()(1);
      verts_ptr->z  = top.getPos()(2);
      verts_ptr->nx = top.getDir()(0);
      verts_ptr->ny = top.getDir()(1);
      verts_ptr->nz = top.getDir()(2);
      verts_ptr++;

      u = M_PI_2 - du;

      su = sin(u);
      cu = cos(u);
      ru = r * cu;

      double v = 0;
      for( int i = 0; i < m2+1; i++, v += dv ) {

        const double sv = sin(v);
        const double cv = cos(v);

        verts_ptr->x  = float(ru*cv);
        verts_ptr->y  = float(ru*sv);
        verts_ptr->z  = float(r*su);
        verts_ptr->nx = float(cu*cv);
        verts_ptr->ny = float(cu*sv);
        verts_ptr->nz = float(su);
        verts_ptr++;
      }


      const Arrow<float,3> bottom = Arrow<float,3>(Point<float,3>(0,0,-float(r)),Vector<float,3>(0,0,-1));
      verts_ptr->x  = bottom.getPos()(0);
      verts_ptr->y  = bottom.getPos()(1);
      verts_ptr->z  = bottom.getPos()(2);
      verts_ptr->nx = bottom.getDir()(0);
      verts_ptr->ny = bottom.getDir()(1);
      verts_ptr->nz = bottom.getDir()(2);
      verts_ptr++;

      u = M_PI_2 - du*(m1-1);
      su = sin(u);
      cu = cos(u);
      ru = r * cu;

      v = 0;
      for( int i = m2; i >= 0; i--, v += dv ) {

        const double sv = sin(v);
        const double cv = cos(v);

        verts_ptr->x  = float(ru*cv);
        verts_ptr->y  = float(ru*sv);
        verts_ptr->z  = float(r*su);
        verts_ptr->nx = float(cu*cv);
        verts_ptr->ny = float(cu*sv);
        verts_ptr->nz = float(su);
        verts_ptr++;
      }
    }


    // Compute body triangle strips on the sphere
    for( int i = 0; i < m1-2; i++ ) {

      const double u1 = M_PI_2 - du*(i+1);
      const double u2 = M_PI_2 - du*(i+2);

      const double su1 = sin(u1);
      const double cu1 = cos(u1);
      const double ru1 = r * cu1;

      const double su2 = sin(u2);
      const double cu2 = cos(u2);
      const double ru2 = r * cu2;

      for( int j = 0; j < m2+1; j++ ) {

        const double v = j * dv;
        const double sv = sin(v);
        const double cv = cos(v);

        verts_ptr->x  = float(ru1*cv);
        verts_ptr->y  = float(ru1*sv);
        verts_ptr->z  = float(r*su1);
        verts_ptr->nx = float(cu1*cv);
        verts_ptr->ny = float(cu1*sv);
        verts_ptr->nz = float(su1);
        verts_ptr++;

        verts_ptr->x  = float(ru2*cv);
        verts_ptr->y  = float(ru2*sv);
        verts_ptr->z  = float(r*su2);
        verts_ptr->nx = float(cu2*cv);
        verts_ptr->ny = float(cu2*sv);
        verts_ptr->nz = float(su2);
        verts_ptr++;
      }
    }

    _vbo.bufferData( verts.getDim() * sizeof(GL::GLVertexNormal), verts.getPtr(), GL_STATIC_DRAW );

  }



} // END namespace GMlib
