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




#ifndef GM_SCENE_VISUALIZERS_PATHTRACKARROWSVISUALIZER_H
#define GM_SCENE_VISUALIZERS_PATHTRACKARROWSVISUALIZER_H


#include "../gmvisualizer.h"

// gmlib
#include <opengl/bufferobjects/gmvertexbufferobject.h>
#include <opengl/bufferobjects/gmindexbufferobject.h>


namespace GMlib {

  class PathTrackArrows;

  class PathTrackArrowsVisualizer : public Visualizer {
    GM_VISUALIZER(PathTrackArrowsVisualizer)
  public:
    PathTrackArrowsVisualizer(Color color = GMcolor::brown());

    void                          render(const SceneObject *obj, const DefaultRenderer *renderer) const override;
    void                          update() override;

    void                          setColor( const Color& col );
    void                          setLineWidth( float line_width = 1.0f );

    void                          reset(std::vector<Point<float,3>>& arrow_path, unsigned int& element, bool& path_full);

  protected:

    std::vector<Point<float,3>>*  _arrow_path;
    unsigned int*                 _element;
    bool*                         _path_full;
    bool                          _full;

    GL::Program                   _prog;
    GL::VertexBufferObject        _vbo;
    Color                         _color;
    float                         _line_width;
    int                           _no_vertices;

  private:
    void                          _fill(GL::GLVertex *ptr);

  }; // END class PathTrackArrowsVisualizer






  //*************************************************************************
  //*************************************************************************
  //******  Include PathTrackArrowsVisualizer class inline functions.  ******
  //******     The other functions are located in the *.cpp file       ******
  //*************************************************************************
  //*************************************************************************



  inline
  void PathTrackArrowsVisualizer::setColor( const Color& col ) {
    _color = col;
  }



  inline
  void PathTrackArrowsVisualizer::setLineWidth(float line_width) {
    _line_width = line_width;
  }



  inline
  void PathTrackArrowsVisualizer::reset(std::vector<Point<float,3>>& arrow_path, unsigned int& element,  bool& path_full) {
    _arrow_path             = &arrow_path;
    _element                = &element;
    _path_full              = &path_full;
    _full                   = false;
  }



  inline
  void PathTrackArrowsVisualizer::_fill(GL::GLVertex *ptr) {

    if(*_path_full)
      for( unsigned int i = (*_element); i < _no_vertices; i++)
         *ptr++ = {static_cast<GLfloat>((*_arrow_path)[i][0]), static_cast<GLfloat>((*_arrow_path)[i][1]), static_cast<GLfloat>((*_arrow_path)[i][2])};
    for( unsigned int i = 0; i < (*_element); i++)
       *ptr++ = {static_cast<GLfloat>((*_arrow_path)[i][0]), static_cast<GLfloat>((*_arrow_path)[i][1]), static_cast<GLfloat>((*_arrow_path)[i][2])};
  }


} // END namespace GMlib



#endif // GM_SCENE_VISUALIZERS_PATHTRACKARROWSVISUALIZER_H
