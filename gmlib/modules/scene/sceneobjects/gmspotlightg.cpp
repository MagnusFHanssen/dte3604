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



#include "gmspotlightg.h"

// stdlib
#include <vector>

// gmlib
#include "../render/gmdefaultrenderer.h"
#include "../../opengl/gmopenglmanager.h"
#include "../../opengl/shaders/gmgeometryshader.h"



namespace GMlib {


  SpotLightG::SpotLightG(): _size(1.0) {

    setMaterial(GMmaterial::polishedSilver());
    init();
    _update(_size);
  }



  SpotLightG::SpotLightG( const Point<float,3>& pos, const Vector<float,3>& dir, Angle cut_off, double size ) : SpotLight( pos, dir, cut_off ), _size(size) {

    setMaterial(GMmaterial::polishedSilver());
    init();
    _update(_size);
  }



  SpotLightG::SpotLightG( const Color& amb, const Color& dif, const Color& spe, const Point<float,3>& pos, const Vector<float,3>& dir, Angle cut_off, double size) : SpotLight( amb, dif, spe, pos, dir, cut_off ), _size(size) {

    setMaterial(GMmaterial::polishedSilver());
     init();
     _update(_size);
  }



  SpotLightG::SpotLightG( const SpotLight& copy ) : SpotLight( copy ), _size(1.0) {

    setMaterial(GMmaterial::polishedSilver());
    init();
    _update(_size);
  }



  SpotLightG::SpotLightG( const SpotLightG& copy) : SpotLight( copy ) {

    _size = copy._size;
    init();
    _update(_size);
  }



  inline
  void SpotLightG::init() {

    _color_prog.acquire("color");
    _shading_prog.acquire("blinn_phong");

    _light_box_geom_vbo.create();
    _light_geom_vbo.create();
  }




  void SpotLightG::_update(double size) {

    _size = size;
    setSurroundingSphere( Sphere<float,3>( Point<float,3>(_size/2, 0.0f, 0.0f ), _size ) );

    // Create geometry
    int sz = 13;
    double ang = M_2PI/(sz-1);
    float a = float(_size*std::sin(_cutoff.getRad()));
    float b = float(_size*std::cos(_cutoff.getRad()));
    std::vector<Point<float,3>>  pt(sz);
    std::vector<Vector<float,3>> no(sz);

    for(int i=0; i<sz; i++){
        pt[i] = { b, a*float(std::sin(i*ang)), a*float(std::cos(i*ang))};
        no[i] = {-a, b*float(std::sin(i*ang)), b*float(std::cos(i*ang))};
        no[i].normalize();
  }

    // the light box geometry
    //-------------------------
    GL::GLVertexNormal light_p;
    DVector<GL::GLVertexNormal> dp_box(2*sz);
    for(int i=0; i<sz; i++) {
      light_p.x   = 0.05*pt[i][0];
      light_p.y   = 0.05*pt[i][1];
      light_p.z   = 0.05*pt[i][2];
      light_p.nx  = no[i][0];
      light_p.ny  = no[i][1];
      light_p.nz  = no[i][2];
      dp_box[2*i] = light_p;

      light_p.x   = pt[i][0];
      light_p.y   = pt[i][1];
      light_p.z   = pt[i][2];
      light_p.nx  = no[i][0];
      light_p.ny  = no[i][1];
      light_p.nz  = no[i][2];
      dp_box[2*i+1] = light_p;
    }
    _light_box_geom_elements = dp_box.getDim();
    _light_box_geom_vbo.bufferData( _light_box_geom_elements * sizeof(GL::GLVertexNormal), dp_box.getPtr(), GL_STATIC_DRAW);

    // the light geometry
    //----------------------
    DVector<GL::GLVertexNormal> dp_light(sz+1);
    light_p.x   = b;
    light_p.y   = 0;
    light_p.z   = 0;
    light_p.nx  = 1.0f;
    light_p.ny  = 0;
    light_p.nz  = 0;
    dp_light[0] = light_p;
    for(int i=0; i<sz; i++) {
      UnitVector<float,3> nor = pt[i];
      light_p.x   = 0.9*pt[i][0];
      light_p.y   = 0.9*pt[i][1];
      light_p.z   = 0.9*pt[i][2];
      light_p.nx  = nor[0];
      light_p.ny  = nor[1];
      light_p.nz  = nor[2];
      dp_light[i+1] = light_p;
    }
    _light_geom_elements = dp_light.getDim();
    _light_geom_vbo.bufferData( _light_geom_elements * sizeof(GL::GLVertexNormal), dp_light.getPtr(), GL_STATIC_DRAW);
  }





  void SpotLightG::setCutOff(const Angle& cut_off) {
    SpotLight::setCutOff(cut_off);
    _update(_size);
  }




  void SpotLightG::localDisplay(const DefaultRenderer* renderer) const {

    const Camera* cam = renderer->getCamera();

    const HqMatrix<float,3> &mvmat = this->getModelViewMatrix(cam);
    const HqMatrix<float,3> &pmat =  this->getProjectionMatrix(cam);

    // Render box surface
    // ---------------------
    _shading_prog.bind(); {
      _shading_prog.uniform( "u_mvmat", mvmat );
      _shading_prog.uniform( "u_mvpmat", pmat * mvmat );

      // Lights
      _shading_prog.bindBufferBase( "DirectionalLights", renderer->getDirectionalLightUBO(), 0);
      _shading_prog.bindBufferBase( "PointLights",       renderer->getPointLightUBO(), 1);
      _shading_prog.bindBufferBase( "SpotLights",        renderer->getSpotLightUBO(), 2);

      // Get Material Data
      const Material &m = this->getMaterial();
      _shading_prog.uniform( "u_mat_amb", m.getAmb() );
      _shading_prog.uniform( "u_mat_dif", m.getDif() );
      _shading_prog.uniform( "u_mat_spc", m.getSpc() );
      _shading_prog.uniform( "u_mat_shi", m.getShininess() );

      GL::AttributeLocation   vert_loc = _shading_prog.getAttributeLocation( "in_vertex" );
      GL::AttributeLocation normal_loc = _shading_prog.getAttributeLocation( "in_normal" );

      _light_box_geom_vbo.bind();
         _light_box_geom_vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE,  sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid*>(0x0));
         _light_box_geom_vbo.enable( normal_loc, 3, GL_FLOAT, GL_TRUE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid*>(sizeof(GL::GLVertex)));
            glDrawArrays( GL_TRIANGLE_STRIP, 0, _light_box_geom_elements );
         _light_box_geom_vbo.disable( normal_loc );
         _light_box_geom_vbo.disable( vert_loc );
      _light_box_geom_vbo.unbind();
    } _shading_prog.unbind();

    // Render light surface
    // ---------------------
    _color_prog.bind(); {
       _color_prog.uniform( "u_mvpmat", pmat * mvmat);
       _color_prog.uniform( "u_color", getAmbient());

       GL::AttributeLocation vert_loc2  = _shading_prog.getAttributeLocation( "in_vertex" );

      _light_geom_vbo.bind();
         _light_geom_vbo.enable( vert_loc2, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid*>(0x0));
            glDrawArrays( GL_TRIANGLE_FAN, 0, _light_geom_elements );
         _light_geom_vbo.disable(vert_loc2);
      _light_geom_vbo.unbind();
    } _color_prog.unbind();

  }




  void SpotLightG::localSelect(const Renderer* renderer, const Color& color) const {

    _color_prog.bind(); {
      // Model view and projection matrices
      _color_prog.uniform("u_mvpmat", this->getModelViewProjectionMatrix(renderer->getCamera()));
      _color_prog.uniform("u_color",  color);

      GL::AttributeLocation vert_loc = _color_prog.getAttributeLocation( "in_vertex" );

      _light_box_geom_vbo.bind();
         _light_box_geom_vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid*>(0x0));
            // Draw top and bottom caps
            glDrawArrays( GL_TRIANGLE_STRIP, 0, _light_box_geom_elements );
         _light_box_geom_vbo.disable(vert_loc);
      _light_box_geom_vbo.unbind();

      _light_geom_vbo.bind();
         _light_geom_vbo.enable( vert_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertexNormal), reinterpret_cast<const GLvoid*>(0x0));
            // Draw top and bottom caps
            glDrawArrays( GL_TRIANGLE_FAN, 0, _light_geom_elements );
         _light_geom_vbo.disable(vert_loc);
      _light_geom_vbo.unbind();
    } _color_prog.unbind();
  }

} // END namespace GMlib

