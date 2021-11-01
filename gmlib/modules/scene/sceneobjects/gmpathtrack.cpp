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




#include "gmpathtrack.h"




namespace GMlib {


  PathTrack::PathTrack( unsigned int max_elements, unsigned int elementstride, const Color& c ) {

    _color          = c;
    _path_full      = false;
    _element_stride = elementstride;
    _element        =_stride_current_element = 0;
    _recent_path.resize( max_elements );
    _viz            = nullptr;

    _local_cs       = false;					// NB!!! Plotting in global coordinates.
    setSurroundingSphere( Sphere<float,3>( Point<float,3>( 0.0f, 0.0f, 0.0f ), 100 ) );
    initVisualizer();
  }




  void PathTrack::localSimulate( double /*dt*/ ) {

    if( ( (_stride_current_element++) % _element_stride ) == 0 ) {

      _recent_path[_element++] = getCenterPos();	//adds the current position to the path.
      if( _element == _recent_path.size() ) {
        _element	= 0;
        _path_full	= true;
      }
      this->setEditDone();
    }
  }




} // END namespace GMlib
