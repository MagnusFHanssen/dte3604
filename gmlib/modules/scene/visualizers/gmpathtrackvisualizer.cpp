/*****************************************************************************
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
*****************************************************************************/


// local
#include "gmpathtrackvisualizer.h"

//GMlib
#include <scene/sceneobjects/gmpathtrack.h>
#include <scene/render/gmdefaultrenderer.h>


namespace GMlib {



PathTrackVisualizer::PathTrackVisualizer(Color color) : _color(color), _line_width(1.0f) {

  _prog.acquire("color");
  _vbo.create();

  _full = false;
}




void PathTrackVisualizer::render(const SceneObject *obj, const DefaultRenderer *renderer) const {

  if(_no_vertices>1) {
    GL_CHECK(::glLineWidth(_line_width));

    _prog.bind(); {
        _prog.uniform( "u_mvpmat", obj->getModelViewProjectionMatrix(renderer->getCamera()));
        _prog.uniform( "u_color", _color) ;
        _prog.uniform( "u_selected", false );

        GL::AttributeLocation vert_loc = _prog.getAttributeLocation( "in_vertex" );

        _vbo.bind();
          _vbo.enable(vert_loc, 3, GL_FLOAT, GL_FALSE, sizeof(GL::GLVertex), reinterpret_cast<const GLvoid*>(0x0));
            glDrawArrays( GL_LINE_STRIP, 0, _no_vertices );
          _vbo.disable(vert_loc);
        _vbo.unbind();
      } _prog.unbind();
  }
}




void PathTrackVisualizer::update() {

    _no_vertices = ((*_path_full) ? static_cast<int>((*_path).size()) : static_cast<int>(*_element));

    const GLsizeiptr data_size = _no_vertices * sizeof(GL::GLVertex);

    if(!_full) _vbo.bufferData( data_size, nullptr, GL_STATIC_DRAW );

    GL::GLVertex *ptr = _vbo.mapBuffer<GL::GLVertex>();
      _fill(ptr);
    _vbo.unmapBuffer();

    if((*_path_full)) _full = true;
}



} // END namespace GMlib
