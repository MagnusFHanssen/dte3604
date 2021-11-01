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



#ifndef GM_SCENE_SCENEOBJECT_PATHTRACKARROWS_H
#define GM_SCENE_SCENEOBJECT_PATHTRACKARROWS_H

#include "../gmsceneobject.h"
#include "../visualizers/gmpathtrackarrowsvisualizer.h"


namespace GMlib {


  /*! \class  PathTrackArrows gmPathTrackArrows.h <gmPathTrackArrows>
   *
   * When this object is attached to another *SceneObject, it will keep track of
   * the path of the object. By setting the element size and the stride between
   * points one can control the behaviour of the PathTrackArrows object.
   */
  class PathTrackArrows : public SceneObject {
    GM_SCENEOBJECT(PathTrackArrows)
  public:
    PathTrackArrows(unsigned int max_elements=200, unsigned int elementstride=10, const Color& c=GMcolor::moccasin());

    void                setColor( const Color& c );
    void                setLineWidth(float line_width);
    void                setArrowLength(float length);
    void                setMaxElement(int max_elements, int elementstride = 10);

  protected:
    void                localSimulate( double dt ) override;
    void                initVisualizer();

  private:

    std::vector<Point<float,3>>   _recent_arrows;

    unsigned int                  _element;
    unsigned int                  _element_stride;
    unsigned int                  _stride_current_element;
    bool                          _path_full;
    Color                         _color;
    float                         _arrow_length;

    bool                          _make;

    PathTrackArrowsVisualizer*    _viz;

  }; // END class





  //*************************************************************************
  //*************************************************************************
  //******     Include PathTrackArrows class inline functions.         ******
  //******     The other functions are located in the *.cpp file       ******
  //*************************************************************************
  //*************************************************************************



  inline
  void PathTrackArrows::setColor( const Color& c ) {
    _color = c;
  }



  inline
  void PathTrackArrows::setLineWidth( float line_width ) {
    if(_viz) _viz->setLineWidth(line_width);
  }



  inline
  void PathTrackArrows::setArrowLength( float length ) {
    _arrow_length = length;
  }



  inline
  void PathTrackArrows::setMaxElement(int max_elements, int elementstride) {
    _path_full      = false;
    _element_stride = elementstride;
    _element        = _stride_current_element = 0;
    _recent_arrows.resize( max_elements );

    if(_viz) _viz->reset(_recent_arrows, _element, _path_full);
  }




  inline
  void PathTrackArrows::initVisualizer() {
      _viz = new PathTrackArrowsVisualizer(_color);
      _viz->reset(_recent_arrows, _element, _path_full);
      SceneObject::insertVisualizer(_viz);
  }

} // END namespace

#endif // GM_SCENE_SCENEOBJECT_PATHTRACKARROWS_H

