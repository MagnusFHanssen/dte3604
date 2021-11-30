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



#include "gmpointlightg.h"

// gmlib
#include "gmsphere3d.h"
#include "../camera/gmcamera.h"
#include "../render/gmdefaultrenderer.h"


namespace GMlib {


PointLightG::PointLightG(float r):_sphere(r) {
    init();
}




PointLightG::PointLightG( const Point<float,3>& pos) : PointLight( pos ), _sphere(0.5f) {
    init();
}




PointLightG::PointLightG(const Color& amb, const Color& dif, const Color& spe, const Point<float,3>& pos) : PointLight(amb,dif,spe,pos), _sphere(0.5f) {
    init();
}




PointLightG::PointLightG( const PointLight& copy ) : PointLight( copy ) {
    init();
}




PointLightG::PointLightG( const PointLightG& copy ) : PointLight( copy ) {
    setSurroundingSphere( copy._sphere );
}




void PointLightG::localDisplay(const DefaultRenderer* renderer) const {
    _sphere.render(this->getModelViewMatrix(renderer->getCamera()), this->getProjectionMatrix(renderer->getCamera()), renderer, this->getMaterial());
}




void PointLightG::localSelect(const Renderer* renderer, const Color& color) const {
    _sphere.render( this->getModelViewProjectionMatrix(renderer->getCamera()), color );
}




void PointLightG::init() {
    this->setMaterial(GMlib::GMmaterial::snow());
    this->setSurroundingSphere( _sphere );
}

} // END namespace GMlib
