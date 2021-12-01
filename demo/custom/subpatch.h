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




#ifndef GM_PARAMETRICS_SURFACES_SUBPATCH_H
#define GM_PARAMETRICS_SURFACES_SUBPATCH_H


#include <vector>

// GMlib includes+
#include <parametrics/gmpsurf.h>


namespace GMlib {


template <typename T>
class SubPatch : public SceneObject {
  GM_SCENEOBJECT(SubPatch)
  public:

  SubPatch( PSurf<T,3>* s, const Point<T,2>& p );
  SubPatch( PSurf<T,3>* s, T u, T v);
  SubPatch( const SubPatch<T>& copy );

  virtual ~SubPatch();

  void         setPrimar(SubPatch<T>* s1, SubPatch<T>* s2);
  void         setSlave()                 { _slave = true;}
  bool         isPrimar()                 { return _primar;}
  void         setNumber(int m)           { _number = m; }
  int          getNumber() const override { return _number; }
  Point<T,2>   getParameter() const       { return _par_point; }
  SubPatch<T>* getSlave(bool first=true)  { return first ? _s1:_s2;}

  const
  HqMatrix<float,3>& getMatr() const      { return _matr;}
  void               setMatr(const HqMatrix<T,3>& m);
//void               setParent(SceneObject* obj);

  void         translate(const Vector<float,3>& trans_vector, bool propagate = true ) override;
  void         translateGlobal(const Vector<float,3>& trans_vector, bool propagate = true ) override;
  void         rotateGlobal(Angle a, const Vector<float,3>& rot_axel, bool propagate = true ) override;
  void         rotateGlobal(Angle a, const Point<float,3>& p,const UnitVector<float,3>& d, bool propagate = true ) override;

  void         translateSlave(const Vector<float,3>& trans_vector);
  void         rotateSlave(Angle a, const Point<float,3>& p,const UnitVector<float,3>& rot_axel);

protected:


  // Protected data for the patch
  PSurf<T,3>*   _s;
  Point<T,2>    _par_point;
  Point<T,3>    _center;
  int           _number;

  bool          _primar;
  bool          _slave;
  SubPatch<T>*  _s1;
  SubPatch<T>*  _s2;

  HqMatrix<T,3> _matr;

  void          primar(const Vector<float,3>& trans_vector);
  void          init();

}; // END class SubPatch






//*****************************************
// Constructors and destructor           **
//*****************************************

template <typename T>
inline
SubPatch<T>::SubPatch( PSurf<T,3>* s, const Point<T,2>& p) {

  _s         = s;
  _par_point = p;
  _number    = -1;
  init();
}



template <typename T>
inline
SubPatch<T>::SubPatch( PSurf<T,3>* s, T u, T v) {

  _s         = s;
  _par_point = {u, v};
  _number    = -1;
  init();
}



template <typename T>
inline
SubPatch<T>::SubPatch( const SubPatch<T>& copy ) {
  _s         = copy._s;
  _par_point = copy._par_point;
  _number    = copy._number;
  init();
}



template <typename T>
SubPatch<T>::~SubPatch() {}




template <typename T>
inline
void SubPatch<T>::setPrimar(SubPatch<T>* s1, SubPatch<T>* s2) {
  _primar = true;
  _s1 = s1;
  _s2 = s2;
  _s1->setSlave();
  _s2->setSlave();
}




template <typename T>
void SubPatch<T>::setMatr( const HqMatrix<T,3>& m) {

  _matr = m;
}




template <typename T>
void SubPatch<T>::translate(const Vector<float,3>& trans_vector, bool propagate ) {

  if(_slave) return;
  SceneObject::translate( trans_vector, propagate );
  _center += trans_vector;
  _matr.translate(trans_vector);
  primar(trans_vector);
}



template <typename T>
void SubPatch<T>::translateGlobal(const Vector<float,3>& trans_vector, bool propagate) {

  if(_slave) return;
  SceneObject::translateGlobal( trans_vector, propagate );
   _matr.translate(trans_vector);
   primar(trans_vector);
}



template <typename T>
void SubPatch<T>::rotateGlobal( Angle a, const Vector<float,3>& rot_axel, bool propagate ) {

  if(_slave) return;
  SceneObject::rotateGlobal( a, rot_axel, propagate );
  _matr = _matrix;
  _matr.translate(-_center);
  if(_primar){
    _s1->setMatr(_matr);
    _s2->setMatr(_matr);
    _s1->rotateSlave(a, getPos(), rot_axel);
    _s2->rotateSlave(a, getPos(), rot_axel);
  }
}



template <typename T>
void SubPatch<T>::rotateGlobal( Angle a, const Point<float,3>& p,const UnitVector<float,3>& d, bool propagate ) {

  if(_slave) return;
  SceneObject::rotateGlobal( a, p, d, propagate );
  _matr = _matrix;
  _matr.translate(-_center);
  if(_primar){
    _s1->setMatr(_matr);
    _s2->setMatr(_matr);
    _s1->rotateSlave(a, getPos()+p, d);
    _s2->rotateSlave(a, getPos()+p, d);
  }
}



template <typename T>
inline
void SubPatch<T>::translateSlave(const Vector<float,3>& trans_vector) {
  SceneObject::translateGlobal( trans_vector, false );
}



template <typename T>
inline
void SubPatch<T>::rotateSlave(Angle a, const Point<float,3>& p,const UnitVector<float,3>& rot_axel) {
  SceneObject::rotateGlobal( a, p, rot_axel, false);
}



template <typename T>
inline
void SubPatch<T>::primar(const Vector<float,3>& trans_vector) {

  if(_primar){
    _s1->setMatr(_matr);
    _s2->setMatr(_matr);
    _s1->translateSlave(trans_vector);
    _s2->translateSlave(trans_vector);
  }
}


template <typename T>
inline
void SubPatch<T>::init() {
  _center    = _s->getPosition(_par_point[0], _par_point[1]);
  SceneObject::translate(_center);
  this->setSurroundingSphere(Sphere<float,3>(Point<T,3>(T(0)), 1.0));

  _primar = false;
  _slave = false;
  _s1 = nullptr;
  _s2 = nullptr;
}

} // END namepace GMlib


#endif // GM_PARAMETRICS_SURFACES_SUBPATCH_H

