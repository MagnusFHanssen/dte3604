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




#ifndef GM_SCENE_VISUALIZERS_SURROUNDINGSPHEREVISUALIZER_H
#define GM_SCENE_VISUALIZERS_SURROUNDINGSPHEREVISUALIZER_H

#include "../gmvisualizer.h"

// gmlib
#include <opengl/bufferobjects/gmvertexbufferobject.h>
#include <opengl/bufferobjects/gmindexbufferobject.h>


namespace GMlib {

  class DefaultRenderer;


  class SurroundingSphereVisualizer : public Visualizer {
      GM_VISUALIZER(SurroundingSphereVisualizer)
    public:
      SurroundingSphereVisualizer(int m1 = 8, int m2 = 8);

      void                      render(const SceneObject *obj, const DefaultRenderer *renderer) const override;
      void                      update() override;

      void                      setSphere(const Sphere<float,3>& s);
      void                      setColor( const Color& color );
      void                      setLineWidth( float line_width = 1.0f );
      void                      setScale(const Vector<float,3>& s);

    protected:
      const int                 _m1;
      const int                 _m2;

      int                       _no_indices;
      float                     _line_width;
      Sphere<float,3>           _sphere;
      bool                      _scaled;
      Vector<float,3>           _scale;

      GL::Program               _prog;
      GL::VertexBufferObject    _vbo;
      GL::IndexBufferObject     _ibo;
      Color                     _color;

  private:

      void              _fillVBO();
      void              _fill(DVector<GL::GLVertex>& dp);

      void              _makeLines();

  }; // END class SurroundingSphereVisualizer




  inline
  void SurroundingSphereVisualizer::setColor(const Color& color) { _color = color; }


} // END namespace GMlib

#endif // GM_SCENE_VISUALIZERS_SURROUNDINGSPHEREVISUALIZER_H
