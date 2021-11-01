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



#ifndef GM_SCENE_SCENEOBJECTS_SPHERE3D_H
#define GM_SCENE_SCENEOBJECTS_SPHERE3D_H


#include "../utils/gmmaterial.h"
//#include "../render/gmdefaultrenderer.h"

#include "../../opengl/gmopengl.h"
#include "../../opengl/gmprogram.h"
#include "../../opengl/bufferobjects/gmvertexbufferobject.h"
#include "../../opengl/bufferobjects/gmindexbufferobject.h"


namespace GMlib {


  /*! \class Sphere3D
   *  \brief Pending Documentation
   *
   *  3D Sphere
   */
  class Sphere3D : public Sphere<float,3> {
  public:
    Sphere3D(float r=1.0, int m1=7, int m2=7);
    Sphere3D(const Sphere<float,3>& s, int m1=7, int m2=7);
    Sphere3D(const Sphere3D& copy );

    void        render(const HqMatrix<float,3>& mvmat, const HqMatrix<float,3>& pmat, const DefaultRenderer *renderer, const Material& mat ) const;

    void        render(const HqMatrix<float,3>& mvpmat, const Color& col ) const;

    void        select(const GL::Program& prog) const;

    void        replot(int m1, int m2);


    void        update(int m1, int m2);


  private:
    GL::Program               _prog;
    GL::Program               _color_prog;

    GL::VertexBufferObject    _vbo;
    GL::IndexBufferObject     _ibo;

    int                       _top_bot_verts;
    int                       _mid_strips;
    int                       _mid_strips_verts;

    int                       _m1, _m2;

    void replot(int m1, int m2, GL::GLVertexNormal * verts_ptr);

  }; // END class Sphere3D








  /*! Sphere3D::Sphere3D(float r, int m1, int m2)
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  Sphere3D::Sphere3D(float r, int m1, int m2) : Sphere<float,3>(Point<float,3>(0.0f),r) {

    _prog.acquire("blinn_phong");
    _color_prog.acquire("color");
    _vbo.create();

    replot(m1,m2);
  }


  /*! Sphere3D::Sphere3D(const Sphere<float,3>& s, int m1, int m2)
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  Sphere3D::Sphere3D(const Sphere<float,3>& s, int m1, int m2) : Sphere<float,3>(s) {

    _prog.acquire("blinn_phong");
    _color_prog.acquire("color");
    _vbo.create();

    replot(m1,m2);
  }

  inline
  Sphere3D::Sphere3D( const Sphere3D& copy ): Sphere<float,3>( copy ) {

    _prog.acquire("blinn_phong");
    _color_prog.acquire("color");
    _vbo.create();

    replot(copy._m1,copy._m2);
  }


  /*! void Sphere3D::display()
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  void Sphere3D::render( const HqMatrix<float,3>& mvmat, const HqMatrix<float,3>& pmat, const DefaultRenderer *renderer, const Material & mat ) const {

    _prog.bind(); {

      // Model view and projection matrices
        _prog.uniform( "u_mvmat", mvmat );
        _prog.uniform( "u_mvpmat", pmat * mvmat );

        // Lights
//        _prog.bindBufferBase( "DirectionalLights",  renderer->getDirectionalLightUBO(), 0 );
//        _prog.bindBufferBase( "PointLights",        renderer->getPointLightUBO(), 1 );
//        _prog.bindBufferBase( "SpotLights",         renderer->getSpotLightUBO(), 2 );

        // Material data
        _prog.uniform( "u_mat_amb", mat.getAmb() );
        _prog.uniform( "u_mat_dif", mat.getDif() );
        _prog.uniform( "u_mat_spc", mat.getSpc() );
        _prog.uniform( "u_mat_shi", mat.getShininess() );

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
    } _color_prog.unbind();
    
  }


  inline
  void Sphere3D::render( const HqMatrix<float,3>& modelview_projection, const Color& color ) const {

    _color_prog.bind(); {

      // Model view and projection matrices
      _color_prog.uniform( "u_mvpmat", modelview_projection );
      _color_prog.uniform( "u_color",  color );

      GL::AttributeLocation vert_loc = _color_prog.getAttributeLocation( "in_vertex" );

      _vbo.bind();
      _vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<const GLvoid*>(0x0) );

      // Draw top and bottom caps
      for( int i = 0; i < 2; i++ )
        glDrawArrays( GL_TRIANGLE_FAN, i * _top_bot_verts, _top_bot_verts );

      // Draw body strips
      for( int i = 0; i < _mid_strips; i++ )
        glDrawArrays( GL_TRIANGLE_STRIP, _top_bot_verts*2 + i*_mid_strips_verts, _mid_strips_verts );

      _vbo.disable(vert_loc);
      _vbo.unbind();

    } _color_prog.unbind();

  }



//  /*! void Sphere3D::display()
//   *  \brief Pending Documentation
//   *
//   *  Pending Documentation
//   */
//  inline
//  void Sphere3D::render( const HqMatrix<float,3>& modelview, const HqMatrix<float,3>& projection, const Material& material  ) const {

//    _shade_prog.bind(); {

//      // Model view and projection matrices
//      _shade_prog.uniform( "u_mvmat",  modelview );
//      _shade_prog.uniform( "u_mvpmat", modelview * projection );

//      // Lights
//      _shade_prog.uniformBlockBinding( "Lights", _lights_ubo, 0 );

//      // Get Material Data
//      _shade_prog.uniform( "u_mat_amb",  material.getAmb() );
//      _shade_prog.uniform( "u_mat_dif",  material.getDif() );
//      _shade_prog.uniform( "u_mat_spc",  material.getSpc() );
//      _shade_prog.uniform( "u_mat_shin", material.getShininess() );

//      GL::AttributeLocation vert_loc = _shade_prog.getAttributeLocation( "in_vertex" );
//      GL::AttributeLocation normal_loc = _shade_prog.getAttributeLocation( "in_normal" );

//      _vbo_v.bind();
//      _vbo_v.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<const GLvoid*>(0x0) );

//      _vbo_n.bind();
//      _vbo_n.enable( normal_loc, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<const GLvoid*>(0x0) );

//      // Draw top and bottom caps
//      for( int i = 0; i < 2; i++ )
//        glDrawArrays( GL_TRIANGLE_FAN, i * _top_bot_verts, _top_bot_verts );

//      // Draw body strips
//      for( int i = 0; i < _mid_strips; i++ )
//        glDrawArrays( GL_TRIANGLE_STRIP, _top_bot_verts*2 + i*_mid_strips_verts, _mid_strips_verts );

//      _vbo_v.bind();
//      _vbo_v.disable(vert_loc);
//      _vbo_n.bind();
//      _vbo_n.disable(normal_loc);
//      _vbo_n.unbind();

//    } _shade_prog.unbind();
//  }

  inline
  void Sphere3D::select( const GL::Program& prog ) const {


    // Model view and projection matrices

    GL::AttributeLocation vert_loc = prog.getAttributeLocation( "in_vertex" );

    _vbo.bind();
    _vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, 0, reinterpret_cast<const GLvoid*>(0x0) );

//    std::cout << "Selthe Sphere3D in color mode; top/bot verts: " << _top_bot_verts << "; mid strips: " << _mid_strips << std::endl;
    // Draw top and bottom caps
    for( int i = 0; i < 2; i++ )
      glDrawArrays( GL_TRIANGLE_FAN, i * _top_bot_verts, _top_bot_verts );

    // Draw body strips
    for( int i = 0; i < _mid_strips; i++ )
      glDrawArrays( GL_TRIANGLE_STRIP, _top_bot_verts*2 + i*_mid_strips_verts, _mid_strips_verts );

    _vbo.disable(vert_loc);
    _vbo.unbind();

  }


  /*! void Sphere3D::replot(int m1, int m2)
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  void Sphere3D::update(int m1, int m2) {

      // Check lower replot boundaries
      if( m1 > 1 )
          _m1 = m1;

      if( m2 > 1 )
          _m2 = m2;

          GL::GLVertexNormal *ptr = _vbo.mapBuffer<GL::GLVertexNormal>();
          replot(m1,  m2, ptr);
          _vbo.unmapBuffer();
}



  /*! void Sphere3D::replot(int m1, int m2)
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  void Sphere3D::replot(int m1, int m2) {

      // Check lower replot boundaries
      if( m1 > 1 )
          _m1 = m1;

      if( m2 > 1 )
          _m2 = m2;

      _top_bot_verts    = _m2+2;
      _mid_strips       = _m1-2;
      _mid_strips_verts = (_m2+1) * 2;

      const unsigned int no_verts = _top_bot_verts * 2 + _mid_strips * _mid_strips_verts;

      DVector<GL::GLVertexNormal> verts(no_verts);
      GL::GLVertexNormal *verts_ptr = verts.getPtr();

      replot(m1,  m2, verts_ptr);

      _vbo.bufferData( verts.getDim() * sizeof(GL::GLVertexNormal), verts.getPtr(), GL_STATIC_DRAW );
}



  /*! void Sphere3D::replot(int m1, int m2, GL::GLVertexNormal * verts_ptr)
   *  \brief Pending Documentation
   *
   *  Pending Documentation
   */
  inline
  void Sphere3D::replot(int m1, int m2, GL::GLVertexNormal * verts_ptr) {

    // Set some vars
    const double x = _pos[0];
    const double y = _pos[1];
    const double z = _pos[2];
    const double r = _value;

    // Compute stride in the spheres u and v parametric direction.
    const double du = M_PI/_m1;
    const double dv = M_2PI/_m2;


    // Compute top triangle fan for the sphere
    const Arrow<float,3> top = Arrow<float,3>(_pos+Point<float,3>(0,0,r), Vector<float,3>(0,0,1));
    verts_ptr->x  = top.getPos()(0);
    verts_ptr->y  = top.getPos()(1);
    verts_ptr->z  = top.getPos()(2);
    verts_ptr->nx = top.getDir()(0);
    verts_ptr->ny = top.getDir()(1);
    verts_ptr->nz = top.getDir()(2);
    verts_ptr++;

    double u  = M_PI_2 - du;
    double v  = 0;
    double su = sin(u);
    double cu = cos(u);
    double ru = r * cu;

    for( int i = 0; i < _m2+1; i++, v += dv ) {
      const double sv = sin(v);
      const double cv = cos(v);

      verts_ptr->x  = float(x + ru*cv);
      verts_ptr->y  = float(y + ru*sv);
      verts_ptr->z  = float(z + r*su);
      verts_ptr->nx = float(cu*cv);
      verts_ptr->ny = float(cu*sv);
      verts_ptr->nz = float(su);
      verts_ptr++;
    }

    // Compute bottom triangle fan for the sphere
    const Arrow<float,3> bottom = Arrow<float,3>(_pos+Point<float,3>(0,0,-r),Vector<float,3>(0,0,-1));
    verts_ptr->x  = bottom.getPos()(0);
    verts_ptr->y  = bottom.getPos()(1);
    verts_ptr->z  = bottom.getPos()(2);
    verts_ptr->nx = bottom.getDir()(0);
    verts_ptr->ny = bottom.getDir()(1);
    verts_ptr->nz = bottom.getDir()(2);
    verts_ptr++;

    u  = M_PI_2 - du*(_m1-1);
    v  = 0;
    su = sin(u);
    cu = cos(u);
    ru = r * cu;

    for( int i = _m2; i >= 0; i--, v += dv ) {
      const double sv = sin(v);
      const double cv = cos(v);

      verts_ptr->x  = float(x + ru*cv);
      verts_ptr->y  = float(y + ru*sv);
      verts_ptr->z  = float(z + r*su);
      verts_ptr->nx = float(cu*cv);
      verts_ptr->ny = float(cu*sv);
      verts_ptr->nz = float(su);
      verts_ptr++;
    }


    // Compute body triangle strips on the sphere
    for( int i = 0; i < _m1-2; i++ ) {

      const double u1 = M_PI_2 - du*(i+1);
      const double u2 = M_PI_2 - du*(i+2);

      const double su1 = sin(u1);
      const double cu1 = cos(u1);
      const double ru1 = r * cu1;

      const double su2 = sin(u2);
      const double cu2 = cos(u2);
      const double ru2 = r * cu2;

      for( int j = 0; j < _m2+1; j++ ) {
        const double v = j * dv;
        const double sv = sin(v);
        const double cv = cos(v);

        verts_ptr->x  = float(x + ru1*cv);
        verts_ptr->y  = float(y + ru1*sv);
        verts_ptr->z  = float(z + r*su1);
        verts_ptr->nx = float(cu1*cv);
        verts_ptr->ny = float(cu1*sv);
        verts_ptr->nz = float(su1);
        verts_ptr++;

        verts_ptr->x  = float(x + ru2*cv);
        verts_ptr->y  = float(y + ru2*sv);
        verts_ptr->z  = float(z + r*su2);
        verts_ptr->nx = float(cu2*cv);
        verts_ptr->ny = float(cu2*sv);
        verts_ptr->nz = float(su2);
        verts_ptr++;
      }
    }
  }


} // END namespace GMlib


#endif // GM_SCENE_SCENEOBJECTS_SPHERE3D_H
