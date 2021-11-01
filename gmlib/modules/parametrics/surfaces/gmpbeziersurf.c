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



#include "../evaluators/gmevaluatorstatic.h"

// gmlib
#include <core/containers/gmdmatrix.h>
#include <scene/selector/gmselector.h>
#include <scene/visualizers/gmselectorgridvisualizer.h>

namespace GMlib {



//*****************************************
// Constructors and destructor           **
//*****************************************

  template <typename T>
  inline
  PBezierSurf<T>::PBezierSurf( const DMatrix< Vector<T, 3> >& cp ) {

    init();

    // Set Control Points
    setControlPoints( cp );
  }


  template <typename T>
  inline
  PBezierSurf<T>::PBezierSurf( const DMatrix< Vector<T, 3> >& c, T s_u, T u, T e_u, T s_v, T v, T e_v  ) {

    init();
    this->setDomainU(s_u, e_u);
    this->setDomainV(s_v, e_v);
    _u_t = (u-s_u)/(e_u-s_u);
    _v_t = (v-s_v)/(e_v-s_v);

    // Generate the control points
    DMatrix<T> bu, bv;
    EvaluatorStatic<T>::evaluateBhp( bu, c.getDim1()-1, ( u - s_u ) / ( e_u - s_u ), T(1)/(e_u-s_u) );
    EvaluatorStatic<T>::evaluateBhp( bv, c.getDim2()-1, ( v - s_v ) / ( e_v - s_v ), T(1)/(e_v-s_v) );
    setScale(T(1)/(e_u-s_u),T(1)/(e_v-s_v));
    bu.invert();
    bv.invert();
    bv.transpose();

    _c = bu * (c^bv);

    for( int i = 0; i < c.getDim1(); i++ )
      for( int j = 0; j < c.getDim2(); j++ )
        _c[i][j] -= c(0)(0);

    this->translateGlobal( c(0)(0) );
  }


  template <typename T>
  inline
  PBezierSurf<T>::PBezierSurf( const PBezierSurf<T>& copy ) : PSurf<T,3>( copy ) {

    init();

    _c = copy._c;
  }


  template <typename T>
  inline
  PBezierSurf<T>::~PBezierSurf() {

    if(_sgv) delete _sgv;
  }


  //*****************************************
  //            Local functons             **
  //*****************************************

  template <typename T>
  inline
  DMatrix< Vector<T,3> >& PBezierSurf<T>::getControlPoints() const {
    return _c;
  }

  template <typename T>
  inline
  int PBezierSurf<T>::getDegreeU() const {
    return _c.getDim1() - 1;
  }


  template <typename T>
  inline
  int PBezierSurf<T>::getDegreeV() const {
    return _c.getDim2() - 1;
  }


  template <typename T>
  inline
  bool PBezierSurf<T>::isSelectorsVisible() const {
    return _selectors;
  }



  template <typename T>
  inline
  void PBezierSurf<T>::setClosed( bool closed_u, bool closed_v ) {

    _cu = closed_u;
    _cv = closed_v;
  }


  template <typename T>
  inline
  void PBezierSurf<T>::setControlPoints( const DMatrix< Vector<T,3> >& cp ) {

    int n1 = _c.getDim1();
    int n2 = _c.getDim2();

    _c = cp;

    if(n1 != _c.getDim1())
      preSample( 1, this->getSamplesU());
    if(n2 != _c.getDim2())
      preSample( 2, this->getSamplesV());

    if(_selectors) {
      hideSelectors();
      showSelectors( _selector_radius, _grid, _selector_color, _grid_color );

    }
  }



  template <typename T>
  inline
  void PBezierSurf<T>::setScale( T su, T sv ) {

      if( su != _su) {
          _su = su;
          preSample( 1, this->getSamplesU());
      }

      if( sv != _sv) {
          _sv = sv;
          preSample( 2, this->getSamplesV());
      }
  }



  template <typename T>
  inline
  void PBezierSurf<T>::updateCoeffs( const Vector<T,3>& d ) {

    if( _c_moved ) {

      HqMatrix<T,3> invmat = this->_matrix;
      invmat.invertOrthoNormal();

      Vector<T,3> diff = invmat*d;
      for( int i = 0; i < _c.getDim1(); i++ ) {
        for( int j = 0; j < _c.getDim2(); j++ ) {

          _c[i][j] += diff;
          _s[i][j]->translateParent( diff );
        }
      }
      this->translateParent( -d, false );
    }
  }




  //********************************************************
  // Overrided (public) virtual functons from SceneObject **
  //********************************************************


  template <typename T>
  void PBezierSurf<T>::edit( int selector_id, const Vector<T,3>& dp  ) {

    _c_moved = true;
      if( this->_parent ) this->_parent->edit( this );
      if( this->_derived ) this->_derived->edit( this );
      _pos_change.push_back(EditSet(mapIndex(selector_id), dp));
      this->setEditDone();
    _c_moved = false;
  }



  // This function is not meant for public use

  /*! void PBezierSurf<T>::replot() const
   *  Not for public use,
   *  For the system to replot after object editing.
   *  We update sampling when control points has been moved.
   *
   */
  template <typename T>
  void PBezierSurf<T>::replot() const {

      updateSamples();
      this->resampleNormals( this->_visu[0][0].sample_val, this->_visu[0][0].normals );
      this->setSurroundingSphere( this->_visu[0][0].sample_val );

      SceneObject::replot();
  }


  //**************************************************
  // Overrided (public) virtual functons from PSurf **
  //**************************************************


  template <typename T>
  bool PBezierSurf<T>::isClosedU() const {
    return _cu;
  }


  template <typename T>
  bool PBezierSurf<T>::isClosedV() const {
    return _cv;
  }



  template <typename T>
  void PBezierSurf<T>::showSelectors( T rad, bool grid, const Color& selector_color, const Color& grid_color ) {

    if( !_selectors ) {
      _s.setDim( _c.getDim1(), _c.getDim2() );

      for( int i = 0, s_id = 0; i < _c.getDim1(); i++ )
        for( int j = 0; j < _c.getDim2(); j++ )
          this->insert( _s[i][j] = new Selector<T,3>( _c[i][j], s_id++, this, rad, selector_color ));
      _selectors = true;
    }
    _selector_radius = rad;
    _selector_color  = selector_color;

    if( grid ) {
      if(!_sgv) _sgv = new SelectorGridVisualizer<T>;
      _sgv->setSelectors( _c );
      _sgv->setColor( grid_color );
      SceneObject::insertVisualizer( _sgv );
      this->setEditDone();
    }
    _grid_color      = grid_color;
    _grid            = grid;
  }



  template <typename T>
  void PBezierSurf<T>::hideSelectors() {

    // Remove Selectors
    if( _selectors ) {
      for( int i = 0; i < _s.getDim1(); i++ ) {
        for( int j = 0; j < _s.getDim2(); j++ ) {
          this->remove( _s[i][j] );
          delete _s[i][j];
        }
      }
      _selectors = false;
    }

    // Hide Selector Grid
    if(_sgv) {
      SceneObject::removeVisualizer( _sgv );
      _sgv->reset();
    }
  }





  /*! void PBezierSurf<T>::toggleSelectors()
   *  To toggle the selectors and selector grid.
   */
  template <typename T>
  void PBezierSurf<T>::toggleSelectors() {

    if(_selectors)  hideSelectors();
    else            showSelectors(_selector_radius, _grid, _selector_color, _grid_color);
  }



  //*****************************************************
  // Overrided (protected) virtual functons from PSurf **
  //*****************************************************

  template <typename T>
  void PBezierSurf<T>::eval( T u, T v, int du, int dv, bool /*lu*/, bool /*lv*/ ) const {
if(u<0 || u>1)
    std::cout << "(u,v) = (" << u << ", " << v << ")" << " U-domene= (" << this->getParStartU() << ", " << this->getParEndU() << ")"<< std::endl;
if(v<0 || v>1)
    std::cout << "(u,v) = (" << u << ", " << v << ")" << " V-domene= (" << this->getParStartV() << ", " << this->getParEndV() << ")"<< std::endl;

      // Set Dimensions
      this->_p.setDim( du+1, dv+1 );

      DMatrix<T> bu, bv;
      EvaluatorStatic<T>::evaluateBhp( bu, this->getDegreeU(), u, _su );
      EvaluatorStatic<T>::evaluateBhp( bv, this->getDegreeV(), v, _sv );

      multEval( bu, bv, du, dv);
  }


  template <typename T>
  T PBezierSurf<T>::getStartPU() const {
    return T(0);
  }


  template <typename T>
  T PBezierSurf<T>::getEndPU() const {
    return T(1);
  }


  template <typename T>
  T PBezierSurf<T>::getStartPV() const {
    return T(0);
  }


  template <typename T>
  T PBezierSurf<T>::getEndPV() const {
    return T(1);
  }



  //*****************************************
  //     Local (protected) functons        **
  //*****************************************

  template <typename T>
  inline
  void PBezierSurf<T>::init() {

    this->_type_id = GM_SO_TYPE_SURFACE_BEZIER;

    _selectors = false;
    _c_moved   = false;

    _cu = false;
    _cv = false;
    _su = T(1);
    _sv = T(1);

    _sgv = 0x0;
  }



  /*! void  PBezierSurf<T>::updateSamples() const
   *  Protected,
   *  Updating sample points and derivatives when control points has been moved,
   */
  template <typename T>
  void  PBezierSurf<T>::updateSamples() const {

      while(_pos_change.size()>0) {
          EditSet es = _pos_change.back();
          for(int i=_pos_change.size()-2; i>=0; i--)
              if (_pos_change[i].ind == es.ind) {
                  es.dp += _pos_change[i].dp;
                  _pos_change.erase(_pos_change.begin()+i);
              }
          for(uint i=0; i<_pre_u.size(); i++)
            for(uint j=0; j<_pre_v.size(); j++)
              comp(this->_visu[0][0].sample_val[i][j], _pre_u[i], _pre_v[j], es.dp, es.ind);
          _pos_change.pop_back();
//          if(_local_pre_eval) {
//              resample( this->_visu[1].sample_val, this->_visu[1].sur_sphere , this->_visu[1], this->_visu[1].sample_val[0].getDim()-1);
//              this->setEditDone();
//          }
      }
      setSurroundingSphere(this->_visu[0][0].sample_val);
  }



  //******************************************
  // Private Overrided  functons from PSurf **
  //******************************************

  template <typename T>
  void PBezierSurf<T>::resample( DMatrix< DMatrix < Vector<T,3> > >& p,
                                 int m1, int m2, int d1, int d2, T s_u, T s_v, T e_u, T e_v ) const{
      // Set Dimensions
      this->_p.setDim(d1+1,d2+1);
      p.setDim(m1, m2);

      for(int i=0; i<m1; i++)
          for(int j=0; j<m2; j++) {
              multEval( _pre_u[i], _pre_v[j], d1, d2);
              p[i][j] = this->_p;
          }
  }


  template <typename T>
  void  PBezierSurf<T>::preSample( int dir, int m ) {
      if( dir==1 )
          internalPreSample( _pre_u, m, _c.getDim1()-1, _su, this->getStartPU(), this->getEndPU() );
      if( dir==2 )
          internalPreSample( _pre_v, m, _c.getDim2()-1, _sv, this->getStartPV(), this->getEndPV() );
  }




  //*****************************
  //  Private Local functons   **
  //*****************************

  template <typename T>
  inline
  void PBezierSurf<T>::internalPreSample( std::vector< DMatrix< T > >& p, int m, int d, T scale, T start, T end ) {

    // compute dt (step in parameter)
    const T dt = ( end - start ) / T(m-1);
    // Set the dimension of the Bernstein-Hermite Polynomial DVector
    p.resize(m);

    // Compute the Bernstein-Hermite Polynomiale, for the B-spline Surface
    for( int j = 0; j < m; j++ )
       EvaluatorStatic<T>::evaluateBhp( p[j], d, j*dt, scale );
  }



  template <typename T>
  inline
  void PBezierSurf<T>::multEval(const DMatrix<T>& bu, const DMatrix<T>& bv, int du, int dv) const {

      int ku = this->getDegreeU()+1;
      int kv = this->getDegreeV()+1;

      DMatrix<Vector<T,3>> c(ku, dv+1);

      // We do these two operations manually here!
      //    bv.transpose();
      //    this->_p = bu * (c^bv);

      //    c= _c^bvT
      for(int i=0; i< ku; i++)
          for(int j=0; j<=dv; j++) {
              c[i][j] = _c(i)(0)*bv(j)(0);
              for(int k=1; k<kv; k++)
                  c[i][j] += _c(i)(k)*bv(j)(k);
          }
      //    _p = bu * c
      for(int i=0; i<=dv; i++)
          for(int j=0; j<=du; j++) {
              this->_p[i][j] = bu(i)(0)*c[0][j];
              for(int k=1; k<ku; k++)
                  this->_p[i][j] += bu(i)(k)*c[k][j];
          }
  }




  /*! void PBezierSurf<T>::comp(DVector<Vector<T,3>>& p, const DMatrix<T>& B, const Vector<T,3>& c, int k) const
   *  Protected,
   *  Partial vector-matrix computation, ie. actually a vector-vector innerproduct.
   *  where vi compute a vector vith a column vector with index k in the matrix B,
   *  i.e. p = v * B_column[k]
   *
   *  \param[out]  p  Updating the position and d derivatives
   *  \param[in]   Bu The Bernstein-Hermite matrix in u - direction.
   *  \param[in]   Bv The Bernstein-Hermite matrix in v - direction.
   *  \param[in]   c  The distance vector, one control point has been moved
   *  \param[in]   ku the u-index of the column to compute with, also the index of the control point that has been moved
   *  \param[in]   ku the v-index of the column to compute with, also the index of the control point that has been moved
   */
  template <typename T>
  inline
  void  PBezierSurf<T>::comp(DMatrix<Vector<T,3>>& p, const DMatrix<T>& Bu, const DMatrix<T>& Bv, const Vector<T,3>& c, Point<int,2> k) const {

      for(int i=0; i<p.getDim1(); i++)
        for(int j=0; j<p.getDim2(); j++)
          p[i][j] += (Bu(i)(k[0])*Bv(j)(k[1]))*c;
  }



  template <typename T>
  inline
  Point<int,2> PBezierSurf<T>::mapIndex(int i) {

    Point<int,2> r;
    r[0] = i/_c.getDim2();
    r[1] = i - r[0]*_c.getDim2();
    return r;
  }



} // END namespace GMlib


