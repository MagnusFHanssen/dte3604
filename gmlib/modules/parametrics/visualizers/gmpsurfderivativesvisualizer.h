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





#ifndef GM_PARAMETRICS_VISUALIZERS_PSURFDERIVATIVESVISUALIZER_H
#define GM_PARAMETRICS_VISUALIZERS_PSURFDERIVATIVESVISUALIZER_H


#include "gmpsurfvisualizer.h"


// gmlib
#include <core/types/gmpoint.h>
#include <core/containers/gmdvector.h>
#include <core/containers/gmdmatrix.h>
#include <core/utils/gmcolor.h>






namespace GMlib {

  enum GM_SURF_DERIVATIVESVISUALIZER_SIZE {
    GM_SURF_DERIVATIVESVISUALIZER_NORMALIZED,
    GM_SURF_DERIVATIVESVISUALIZER_RELATIVE,
    GM_SURF_DERIVATIVESVISUALIZER_ABSOLUTE
  };

  template <typename T, int n>
  class PSurfDerivativesVisualizer : public PSurfVisualizer<T,n> {
    GM_VISUALIZER(PSurfDerivativesVisualizer)
  public:
    PSurfDerivativesVisualizer(int u=1, int v=0);
    PSurfDerivativesVisualizer( const PSurfDerivativesVisualizer<T,n>& copy );

    virtual ~PSurfDerivativesVisualizer();

    void          render(const SceneObject *obj, const DefaultRenderer *renderer) const;
    const Color&  getColor() const;
    int           getDerivativeU() const;
    int           getDerivativeV() const;

    GM_SURF_DERIVATIVESVISUALIZER_SIZE getMode() const;

    double        getSize() const;

    void          update() override;

    void          setColor( const Color& color );
    void          setDerivatives( int u = 0, int v = 0 );
    void          setMode( GM_SURF_DERIVATIVESVISUALIZER_SIZE mode );
    void          setSize( double size = 1.0 );

  protected:
    Color                               _color;
    GL::Program                         _prog;
    GL::VertexBufferObject              _vbo;
    int                                 _no_elements;
    GM_SURF_DERIVATIVESVISUALIZER_SIZE  _mode;
    double                              _size;
    int                                 _u;
    int                                 _v;
  };


} // END namespace GMlib

// Include PSurfDerivativesVisualizer class function implementations
#include "gmpsurfderivativesvisualizer.c"

#endif // GM_PARAMETRICS_VISUALIZERS_PSURFDERIVATIVESVISUALIZER_H
