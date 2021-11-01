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





// gmlib
#include "visualizers/gmpsurfdefaultvisualizer.h"
#include <core/utils/gmdivideddifferences.h>


// stl
#include <cmath>
#include <sstream>
#include <iomanip>

namespace GMlib {

    //*************************************
    //***  Constructors and destructors  **
    //*************************************

  template <typename T, int n>
  inline
  PSurf<T,n>::PSurf( int s1, int s2 ) {

    _no_sam_u                       = s1;
    _no_sam_v                       = s2;
    _no_der_u                       = 1;
    _no_der_v                       = 1;

    _d1                             = -1;
    _d2                             = -1;
    _tr_u                           = T(0);
    _sc_u                           = T(1);
    _is_scaled_u                    = false;
    _tr_v                           = T(0);
    _sc_v                           = T(1);
    _is_scaled_v                    = false;

    _resample                       = false;

    setNoDer( 2 );

//    _default_visualizer = 0x0;
  }


  template <typename T, int n>
  inline
  PSurf<T,n>::PSurf( const PSurf<T,n>& copy ) : Parametrics<T,2,n>( copy ) {

    _p            = copy._p;
    _u            = copy._u;
    _v            = copy._v;
    _d1           = copy._d1;
    _d2           = copy._d2;
    _tr_u         = copy._tr_u;
    _sc_u         = copy._sc_u;
    _is_scaled_u  = copy._is_scaled_u;
    _tr_v         = copy._tr_v;
    _sc_v         = copy._sc_v;
    _is_scaled_v  = copy._is_scaled_v;
    _sam_p_u      = copy._sam_p_u;
    _sam_p_v      = copy._sam_p_v;
    _no_sam_p_u   = copy._no_sam_p_u;
    _no_sam_p_v   = copy._no_sam_p_v;
    _default_d    = copy._default_d;

    _no_sam_u     = copy._no_sam_u;
    _no_sam_v     = copy._no_sam_v;
    _no_der_u     = copy._no_sam_u;
    _no_der_v     = copy._no_sam_v;

    _resample     = false;

//    _default_visualizer = 0x0;
  }

  template <typename T, int n>
  PSurf<T,n>::~PSurf() {

    enableDefaultVisualizer( false );
    if( _visu.default_visualizer )
      delete _visu.default_visualizer;
  }







  //******************************************************
  //      public evaluate functons for PSurf           **
  //******************************************************



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluate( T u, T v, int d1, int d2 ) const {

    _eval(u, v, d1, d2);
    return _p;
  }



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluateParent( T u, T v, int d1, int d2 ) const {

    _eval(u,v,d1,d2);
    p.setDim( _p.getDim1(), _p.getDim2() );
    _mat = this->_matrix.template toType<T>();

    p[0][0] = _mat * _p[0][0].toPoint();
    for( int j = 1; j < p.getDim2(); j++ )
      p[0][j] = _mat * _p[0][j];
    for( int i = 1; i < p.getDim1(); i++ )
      for( int j = 0; j < p.getDim2(); j++ )
        p[i][j] = _mat * _p[i][j];

    return p;
  }




  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluateGlobal( T u, T v, int d1, int d2 ) const {

    _eval(u,v,d1,d2);
    p.setDim( _p.getDim1(), _p.getDim2() );
    _mat = this->_present.template toType<T>();

    p[0][0] = _mat * _p[0][0].toPoint();
    for( int j = 1; j < p.getDim2(); j++ )
      p[0][j] = _mat * _p[0][j];
    for( int i = 1; i < p.getDim1(); i++ )
      for( int j = 0; j < p.getDim2(); j++ )
        p[i][j] = _mat * _p[i][j];

    return p;
  }




  template <typename T, int n>
  inline
  DVector<Vector<T,n>>& PSurf<T,n>::evaluateD( T u, T v, int d ) const {
    // Here we are copy  the matrix diagonally into a vector
    static DVector<Vector<T,n> > p;

    _eval(u, v, d, d);
    p.setDim((d*d+3*d+2)/2);

    for(int i = 0, k=0; i <= d; i++)
      for(int j = 0; j <= i; j++)
        p[k++] = _p[i-j][j];

    return p;
  }




  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluate( const Point<T,2>& p, const Point<int,2>& d ) const {

    _eval( p(0), p(1), d(0), d(1) );
    return _p;
  }



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluateParent( const Point<T,2>& p, const Point<int,2>& d ) const {

    return evaluateParent( p(0), p(1), d(0), d(1) );
  }



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluateGlobal( const Point<T,2>& p, const Point<int,2>& d ) const {

    return evaluateGlobal( p(0), p(1), d(0), d(1) );
  }


  template <typename T, int n>
  inline
  DVector<Vector<T,n>>& PSurf<T,n>::evaluateD( const Point<T,2>& p, int d ) const {

    return evaluateD( p(0), p(1), d );
  }



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluate( int i, int j ) const {

    return _pre_val[i][j];
  }



  template <typename T, int n>
  inline
  DMatrix<Vector<T,n>>& PSurf<T,n>::evaluateParent(  int i, int j  ) const {

    static DMatrix<Vector<T,n> > p;
    DMatrix<Vector<T,n>>& q = _pre_val[i][j];
    int k1 = q.getDim1();
    int k2 = q.getDim2();
    p.setDim(k1, k2);

    p[0][0] = _mat * q[0][0].toPoint();
    for( int j = 1; j < k2; j++ )
      p[0][j] = _mat * q[0][j];
    for( int i = 1; i < k1; i++ )
      for( int j = 0; j < k2; j++ )
        p[i][j] = _mat * q[i][j];

    return p;
  }




  //******************************************************
  //      public closest point functions                **
  //******************************************************


  template <typename T, int n>
  inline
  void PSurf<T,n>::estimateClpPar( const Point<T,n>& p, T& u, T& v, int m ) const {

    T su = getParStartU();
    T sv = getParStartV();
    T du = getParDeltaU()/(m-1);
    T dv = getParDeltaV()/(m-1);

    Vector<T,n> q = getPosition(su, sv);
    T min = (p-q).getLength();
    u = su; v = sv;

    for(int i=0; i<m; i++) {
      for(int j=0; j<m; j++) {
        if(!(i==0 && j==0)) {
          q = getPosition(su+i*du, sv+j*dv);
          T mn = (p-q).getLength();
          if (mn < min) {
            min = mn;
            u = su + i*du;
            v = sv + j*dv;
          }
        }
      }
    }
  }




  template <typename T, int n>
  bool PSurf<T,n>::getClosestPoint( const Point<T,n>& q, T& u, T& v, double eps, int max_iterations ) const {

    T a11, a12, a21, a22, b1, b2;
    T du, dv, det;

    HqMatrix<T,n> invmat = this->_present.template toType<T>();
    invmat.invertOrthoNormal();
    Point<T,n> p = invmat * q;

    for(int i = 0; i < max_iterations; i++ ) {

      DMatrix< Vector<T,n>>& r = evaluate( u, v, 2, 2 );
      Vector<T,n> d = p-r[0][0];

      a11 =       d*r[2][0] - r[1][0] * r[1][0];
      a21 = a12 = d*r[1][1] - r[1][0] * r[0][1];
      a22 =       d*r[0][2] - r[0][1] * r[0][1];

      b1  = -(d*r[1][0]);
      b2  = -(d*r[0][1]);
      det = a11*a22 - a12*a21;

      du  = (b1*a22 - a12*b2) / det;
      dv  = (a11*b2 - b1*a21) / det;
      u   += du;
      v   += dv;

      if(std::abs(du) < eps && std::abs(dv) < eps)
        return true;
    }
    return false;
  }



  template <typename T, int n>
  inline
  bool PSurf<T,n>::getClosestPoint( const Point<T,n>& q, Point<T,2>& uv, double eps, int max_iterations ) const {
    return getClosestPoint(q, uv[0], uv[1], eps, max_iterations);
  }





  //**************************************************
  //      public curvature functions                **
  //**************************************************


  template <typename T, int n>
  T PSurf<T,n>::getCurvatureGauss( T u, T v ) const {

    T E, F, G, e, f, g;
    _computeEFGefg( u, v, E, F, G, e, f, g );

    return (e*g - f*f) / (E*G - F*F);
  }


  template <typename T, int n>
  T PSurf<T,n>::getCurvatureMean( T u, T v ) const {

      T E, F, G, e, f, g;
      _computeEFGefg( u, v, E, F, G, e, f, g );

    return 0.5 * (e*G - 2 * (f*F) + g*E) / (E*G - F*F);
  }


  template <typename T, int n>
  T PSurf<T,n>::getCurvaturePrincipalMax( T u, T v ) const {

    T K = getCurvatureGauss( u, v );
    T H = getCurvatureMean( u, v );

    return H + sqrt( H*H - K );
  }


  template <typename T, int n>
  T PSurf<T,n>::getCurvaturePrincipalMin( T u, T v ) const {

    T K = getCurvatureGauss( u, v );
    T H = getCurvatureMean( u, v );

    return H - sqrt( H*H - K );
  }



  //***********************************************************
  //   To see the number of derivatives in pre-evaluation
  //***********************************************************


  template <typename T, int n>
  inline
  int PSurf<T,n>::getDerivativesU() const {

    return _no_der_u;
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getDerivativesV() const {

    return _no_der_v;
  }



  //****************************************************************
  //   To get the position or given derivatives at a given (u,v)
  //****************************************************************


  template <typename T, int n>
  inline
  const Point<T,n>& PSurf<T,n>::operator () ( T u, T v ) const {

    _eval(u, v, _default_d, _default_d);
    return _p[0][0].toPoint();
  }


  template <typename T, int n>
  inline
  const Point<T,n>& PSurf<T,n>::getPosition( T u, T v ) const {

    _eval(u, v, 0, 0);
    return _p[0][0];
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getDerU( T u, T v ) const {

    _eval(u, v, 1, 0);
    return _p(1)(0);
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getDerUU( T u, T v ) const {

    _eval(u, v, 2, 0);
    return _p(2)(0);
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getDerUV( T u, T v ) const {

    _eval(u, v, 2, 2);
    return _p(1)(1);
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getDerV( T u, T v ) const {

    _eval(u, v, 0, 1);
    return _p(0)(1);
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getDerVV( T u, T v ) const {

    _eval(u, v, 0, 2);
    return _p(0)(2);
  }


  template <typename T, int n>
  inline
  const Vector<T,n>& PSurf<T,n>::getNormal() const {

    return _n = _p(1)(0)^_p(0)(1);
  }



  //***********************************************************
  //   To get the desired domain of the surface
  //***********************************************************


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParStartU() const {

    return getStartPU() + _tr_u;
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParEndU() const {

    return getParStartU() + getParDeltaU();
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParDeltaU() const {

    return _sc_u * (getEndPU() - getStartPU());
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParStartV() const {

    return getStartPV() + _tr_v;
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParEndV() const {

    return getParStartV() + getParDeltaV();
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::getParDeltaV() const {

    return _sc_v * (getEndPV() - getStartPV());
  }




  //***********************************************************
  //   To get the sample data
  //***********************************************************


  template <typename T, int n>
  inline
  int PSurf<T,n>::getNumSamIntPU() const {

    return _no_sam_p_u.getDim();
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getNumSamIntPV() const {

    return _no_sam_p_v.getDim();
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getSamPU( int i ) const {

    return _no_sam_p_u(i);
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getSamPV( int i ) const {

    return _no_sam_p_v(i);
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getSamplesU() const {

    return _no_sam_u;
  }


  template <typename T, int n>
  inline
  int PSurf<T,n>::getSamplesV() const {

    return _no_sam_v;
  }



  //******************************************
  // Virtual functons for surface properies **
  //******************************************

  template <typename T, int n>
  bool PSurf<T,n>::isClosedU() const {
    return false;
  }


  template <typename T, int n>
  bool PSurf<T,n>::isClosedV() const {
    return false;
  }



  //*******************************************************
  //***  Virtual functons for pre-samling  and plotting  **
  //*******************************************************

  template <typename T, int n>
  void PSurf<T,n>::preSample( int /*m1*/, int /*m2*/, int /*d1*/, int /*d2*/,
                                T /*s_u*/, T /*s_v*/, T /*e_u*/, T /*e_v*/ ) {}


  template <typename T, int n>
  void PSurf<T,n>::preSample( int /*dir*/, int /*m*/ ) {}



  /*! void  PSurf<T,n>::sample( int m, int d )
   *  This function is the sample function for static surfaces plotted with one partition
   *
   *  \param[in] m1 The number of samples in u - direction
   *  \param[in] m2 The number of samples in v - direction
   *  \param[in] d1 The number of derivatives at each sample in u - direction
   *  \param[in] d2 The number of derivatives at each sample in v - direction
   */
  template <typename T, int n>
  void PSurf<T,n>::sample( int m1, int m2, int d1, int d2 ) {

    initSample(m1, m2, d1, d2);
    _visu[0][0] = {m1,m2};
    _visu[0][0].s_e_u = { getStartPU(), getEndPU()};
    _visu[0][0].s_e_v = { getStartPV(), getEndPV()};
    DMatrix<DMatrix<Vector<T,n>>>&  p = _visu[0][0].sample_val;
    DMatrix<Vector<float,n>>& normals = _visu[0][0].normals;
    Sphere<T,3>&                    s = _visu[0][0].sur_sphere;
    // Calculate sample positions, derivatives, normals and surrounding sphere
    resample( p, m1, m2, d1, d2, _visu[0][0].s_e_u[0], _visu[0][0].s_e_v[0], _visu[0][0].s_e_u[1], _visu[0][0].s_e_v[1] );
    resampleNormals( p, normals );
    computeSurroundingSphere(p, s);

    // Replot Visaulizers
    for( int i = 0; i < _visu[0][0].vis.size(); i++ ) {
      _visu[0][0].vis[i]->set( p );
      _visu[0][0].vis[i]->set( normals );
      _visu[0][0].vis[i]->set( isClosedU(), isClosedV() );
    }
    this->setEditDone();
  }




  template <typename T, int n>
  void PSurf<T,n>::initSample( int& mu, int& mv, int& du, int& dv )  {

      if( mu != _visu.no_sample[0] && mu > 1) { // u - direction
          _visu.no_sample[0] = mu;
          preSample(1, mu);
      }
      else mu = _visu.no_sample[0];

      if( mv != _visu.no_sample[1] && mv > 1) { // v - direction
          _visu.no_sample[1] = mv;
          preSample(2, mv);
      }
      else mv = _visu.no_sample[1];

      // Correct derivatives ?
      if( du < 1 )    du = _visu.no_derivatives[0];
      else            _visu.no_derivatives[0] = du;
      if( dv < 1 )    dv = _visu.no_derivatives[1];
      else            _visu.no_derivatives[1] = dv;
  }




  template <typename T, int n>
  void PSurf<T,n>::replot() const {

    // Give reference to updated data to all visualizers
    for( int i=0; i<_visu.getDim1(); i++ )
      for( int j=0; j<_visu.getDim2(); j++ )
        for(unsigned int r=0; r<_visu[i][j].vis.size(); r++) {
          _visu[i][j].vis[r]->set(_visu[i][j].sample_val);
          _visu[i][j].vis[r]->set(_visu[i][j].normals);
          _visu[i][j].vis[r]->set(_visu.getDim1()>1? false:isClosedU(), _visu.getDim2()>1? false:isClosedV());
        }
    uppdateSurroundingSphere();
    SceneObject::replot();
  }





  //*******************************************************
  //***  functons to handle visualizers to the surface   **
  //*******************************************************

  template <typename T, int n>
  void PSurf<T,n>::enableDefaultVisualizer( bool enable ) {

    if( !enable )
      removeVisualizer(  _visu.default_visualizer );
    else {
      if( !_visu.default_visualizer ){
        _visu.default_visualizer = new PSurfDefaultVisualizer<T,n>(_visu[0][0].sample_val, _visu[0][0].normals);
        _visu.default_visualizer->set(isClosedU(), isClosedV());
      }
      insertVisualizer( _visu.default_visualizer );
    }
  }




  template <typename T, int n>
  PSurfVisualizer<T, n> *PSurf<T,n>::getDefaultVisualizer() const {

    return _visu.default_visualizer;
  }



  template <typename T, int n>
  void PSurf<T,n>::toggleDefaultVisualizer() {

    if( _visu.default_visualizer == nullptr)
      enableDefaultVisualizer( true );
    else
      enableDefaultVisualizer( false );
  }



  template <typename T, int n>
  inline
  void PSurf<T,n>::insertVisualizer( Visualizer *visualizer ) {

    // Is it already active
    if(this->_visualizers.exist( visualizer)) return;

    SceneObject::insertVisualizer( visualizer );

    PSurfVisualizer<T,n> *visu = dynamic_cast<PSurfVisualizer<T,n>*>( visualizer );
    if( !visu ) {
      SceneObject::insertVisualizer( visualizer );
      return;
    }

    // Is this a default visualiser that is not enabled
    if(_visu.default_visualizer == nullptr){
      PSurfDefaultVisualizer<T,n>* visu = dynamic_cast<PSurfDefaultVisualizer<T,n>*>(visualizer);
      if(visu)
        _visu.default_visualizer = visu;
    }

    // If it is a surface visualiser we need to know if it only must be made active or made in total
    int k = -1;
    for(uint i=0; i<_visu[0][0].vis.size(); i++)
      if(visu == _visu[0][0].vis[i])   k=i;

    if(k>=0){
        for(uint i=0; i<_visu.getDim1(); i++)
            for(uint j=0; j<_visu.getDim2(); j++)
            SceneObject::insertVisualizer(_visu[i][j].vis[k]);
    } else{
        _visu[0][0].vis.push_back(visu);
        SceneObject::insertVisualizer(visu);
        for(uint i=1; i<_visu.getDim1(); i++)
            for(uint j=0; j<_visu.getDim2(); j++) {
                PSurfVisualizer<T,n>* v = dynamic_cast<PSurfVisualizer<T,n>*>(visu->makeCopy());
                _visu[i][j].vis.push_back(v);
                SceneObject::insertVisualizer(v);
            }
    }
  }



  template <typename T, int n>
  inline
  void PSurf<T,n>::removeVisualizer( Visualizer *visualizer ) {

    SceneObject::removeVisualizer( visualizer );
  }




  template <typename T, int n>
  void PSurf<T,n>::resample( DMatrix< DMatrix < Vector<T,n> > >& p,
                                    int m1, int m2, int d1, int d2, T s_u, T s_v, T e_u, T e_v ) const {
    _resample = true;
    p.setDim(m1, m2);

    T du = (e_u-s_u)/(m1-1);
    T dv = (e_v-s_v)/(m2-1);

    for(int i=0; i<m1-1; i++) {
      _ind[0]=i;
      T u = s_u + i*du;
      for(int j=0;j<m2-1;j++) {
        _ind[1]=j;
        eval(u, s_v + j*dv, d1, d2, true, true );
        p[i][j] = _p;
      }
      _ind[1]=m2-1;
      eval(u, e_v, d1, d2, true, false);
      p[i][m2-1] = _p;
    }
    _ind[0]=m1-1;
    for(int j=0;j<m2-1;j++) {
      _ind[1]=j;
      eval(e_u, s_v + j*dv, d1, d2, false, true);
      p[m1-1][j] = _p;
    }
    _ind[1]=m2-1;
    eval(e_u, e_v, d1, d2, false, false);
    p[m1-1][m2-1] = _p;

    switch( this->_dm ) {
      case GM_DERIVATION_EXPLICIT:
        // Do nothing, evaluator algorithms for explicite calculation of derivatives
        // should be defined in the eval( ... ) function enclosed by
        // if( this->_derivation_method == this->EXPLICIT ) { ... eval algorithms for derivatives ... }
        break;
      case GM_DERIVATION_DD:
        DD::compute2D(p,double(du),double(dv),isClosedU(),isClosedV(),d1,d2);
        break;
    }

    _resample = false;
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::resample( DMatrix<DMatrix <DMatrix <Vector<T,n> > > >& a,
                                                int m1, int m2, int d1, int d2 ) {

    resample( a, m1, m2, d1, d2, getStartPU(), getStartPV(), getEndPU(), getEndPV() );
  }


  template <typename T, int n>
  void PSurf<T,n>::resampleNormals( const DMatrix<DMatrix<Vector<T,n> > > &p, DMatrix<Vector<float,3> > &normals ) const {

    normals.setDim( p.getDim1(), p.getDim2() );

    for( int i = 0; i < p.getDim1(); i++ )
      for( int j = 0; j < p.getDim2(); j++ ){
        normals[i][j] = p(i)(j)(1)(0) ^ p(i)(j)(0)(1);
        normals[i][j].normalize();
      }
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainU( T start, T end ) {

      setDomainUScale((end - start)/(getEndPU()-getStartPU()));
      setDomainUTrans(start - getStartPU());
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainV( T start, T end ) {

      setDomainVScale((end - start)/(getEndPV()-getStartPV()));
      setDomainVTrans(start - getStartPV());
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainUScale( T sc ) {

    _sc_u = sc;
    if(GMutils::compValueF(sc, T(1)))
      _is_scaled_u = false;
    else
      _is_scaled_u = true;
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainVScale( T sc ) {

    _sc_v = sc;
    if(GMutils::compValueF(sc, T(1)))
      _is_scaled_v = false;
    else
      _is_scaled_v = true;
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainUTrans( T tr ) {

    _tr_u = tr;
  }


  template <typename T, int n>
  inline
  void PSurf<T,n>::setDomainVTrans( T tr ) {

    _tr_v = tr;
  }



  template <typename T, int n>
  inline
  void PSurf<T,n>::setNoDer( int d ) {

     _default_d  = d;
  }



  template <typename T, int n>
  void PSurf<T,n>::setSurroundingSphere( const DMatrix< DMatrix< Vector<T,n> > >& p ) const {
    Sphere<T,n>&  s = _visu[0][0].sur_sphere;
    computeSurroundingSphere(p, s);
    SceneObject::setSurroundingSphere(s);
  }



  template <typename T, int n>
  void PSurf<T,n>::computeSurroundingSphere( const DMatrix<DMatrix<Vector<T,n>>>& p, Sphere<T,n>& s ) const {

      int n1 = p.getDim1()-1;
      int m1 = p.getDim2()-1;
      int n2 = n1/2;
      int m2 = m1/2;
      s.reset();
      // center
      s += p(n2)(m2)(0)(0).toPoint();
      // corner
      s += p( 0)( 0)(0)(0).toPoint();
      s += p(n1)(m1)(0)(0).toPoint();
      s += p(n1)( 0)(0)(0).toPoint();
      s += p( 0)(m1)(0)(0).toPoint();
      // midpoints edge
      s += p(n1)(m2)(0)(0).toPoint();
      s += p( 0)(m2)(0)(0).toPoint();
      s += p(n2)(m1)(0)(0).toPoint();
      s += p(n2)( 0)(0)(0).toPoint();
      // quater points
      if(n1>4 && n2>4) {
          int n4 = n1/4;
          int n3 = 3*n4;
          int m4 = m1/4;
          int m3 = 3*m4;
          s += p(n4)(m4)(0)(0).toPoint();
          s += p(n3)(m3)(0)(0).toPoint();
          s += p(n4)(m3)(0)(0).toPoint();
          s += p(n3)(m4)(0)(0).toPoint();
      }
      if(n1>8 && m1>8) {
          n1 /= 8;
          m1 /= 8;
          int n3 = 3*n1;
          int m3 = 3*m1;
          int n5 = 5*n1;
          int m5 = 5*m1;
          int n7 = 7*n1;
          int m7 = 7*m1;
          s += p(n1)(m1)(0)(0).toPoint();
          s += p(n7)(m7)(0)(0).toPoint();
          s += p(n1)(m7)(0)(0).toPoint();
          s += p(n7)(m1)(0)(0).toPoint();

          s += p(n3)(m3)(0)(0).toPoint();
          s += p(n5)(m5)(0)(0).toPoint();
          s += p(n3)(m5)(0)(0).toPoint();
          s += p(n5)(m3)(0)(0).toPoint();

          s += p(n1)(m3)(0)(0).toPoint();
          s += p(n7)(m5)(0)(0).toPoint();
          s += p(n5)(m7)(0)(0).toPoint();
          s += p(n3)(m1)(0)(0).toPoint();

          s += p(n7)(m3)(0)(0).toPoint();
          s += p(n5)(m1)(0)(0).toPoint();
          s += p(n3)(m7)(0)(0).toPoint();
          s += p(n1)(m5)(0)(0).toPoint();
      }
      if(n1>2 && m1>2)
          for(int i=1; i<p.getDim1()-2; i+=4)
              for(int j=1; j<p.getDim2()-2; j+=4)
                  s += p(i)(j)(0)(0).toPoint();
  }




  /*! void PSurf<T,n>::uppdateSurroundingSphere() const
   *  Summing up all computed surrounding sphere from all partitions
   *  and uppdating the offisial SceenObject surrounding sphere.
   *  NB! This function do not compute the spheres for each partitions.
   */
  template <typename T, int n>
  void PSurf<T,n>::uppdateSurroundingSphere() const {

    // Updating surrounding sphere
    Sphere<T,3> sph;
    for( int i=0; i<_visu.getDim1(); i++ )
      for( int j=0; j<_visu.getDim2(); j++ )
        sph += _visu[i][j].sur_sphere;
    SceneObject::setSurroundingSphere(sph);
  }




  /*! void  PSurf<T,n>::prepareVisualizers()
   *  Private, not for public use
   *  Remove all old visualizers
   *  Creates and inserts new default visualizers for all partitions
   */
  template <typename T, int n>
  void  PSurf<T,n>::prepareVisualizers() {
    // Updating visualizers in new partitions
    for(unsigned int i=0; i < _visu.getDim1(); i++)
      for(unsigned int j=0; j < _visu.getDim2(); j++)
        if(!(i == 0 && j == 0)) {
          if(_visu[i][j].vis.size() == 0) {
            for(uint k=0; k < _visu[0][0].vis.size(); k++) {
              PSurfVisualizer<T,3>* vis = dynamic_cast<PSurfVisualizer<T,3>*>(_visu[0][0].vis[k]->makeCopy());
              _visu[i][j].vis.push_back(vis);
              SceneObject::insertVisualizer(vis);
            }
          }
        }
  }




  /*! void  PSurf<T,n>::cleanVisualizers()
   *  Private, not for public use
   *  Remove all old visualizers exept the one in _visu[0]
   *
   *  \param[in] i    the visualizer to remove and delete.
   *                  If k=0, then all surface-visualisers are deleted
   *                  if k=1, then all exept _visu[0][0] (the original inserted ones) are deleted
   */
  template <typename T, int n>
  void  PSurf<T,n>::cleanVisualizers(int k) {

    // Updating visualizers in new partitions
    for(int i=0; i < _visu.getDim1(); i++)
      for(int j=0; j < _visu.getDim2(); j++)
        if(!(k==1 && i==0 && j==0)) {
          while(_visu[i][j].vis.size()>0) {
            PSurfVisualizer<T,n>* tmp = _visu[i][j].vis.back();
            _visu[i][j].vis.pop_back();
            SceneObject::removeVisualizer(tmp);
            delete tmp;
          }
        }
  }


  template <typename T, int n>
  inline
  T PSurf<T,n>::_mapU( T u ) const {

    return getStartPU() +  ( u - _tr_u - getStartPU() )/_sc_u;
  }



  template <typename T, int n>
  inline
  T PSurf<T,n>::_mapV( T v ) const {

    return getStartPV() + ( v - _tr_v - getStartPV() )/_sc_v;
  }



  template <typename T, int n>
  Parametrics<T,2,n>* PSurf<T,n>::split( T /*t*/, int /*uv*/ ) {

    return 0;
  }



  template <typename T, int n>
  inline
  void PSurf<T,n>::_eval( T u, T v, int d1, int d2 ) const {

    if( !(d1 <= _d1 and d2 <=_d2 and GMutils::compValueF(u,_u) and GMutils::compValueF(v,_v) ) ) {
      _u  = u;
      _v  = v;
      _d1 = d1;
      _d2 = d2;
      eval( _mapU(u), _mapV(v), d1, d2 );
    }
  }



  template <typename T, int n>
  inline
  void PSurf<T,n>::_computeEFGefg( T u, T v, T& E, T& F, T& G, T& e, T& f, T& g ) const {
      _eval(u,v,2,2);
      UnitVector<T,n>  N   = _p[1][0]^_p[0][1];
      Vector<T,n>      du  = _p[1][0];
      Vector<T,n>      dv  = _p[0][1];
      Vector<T,n>      duu = _p[2][0];
      Vector<T,n>      duv = _p[1][1];
      Vector<T,n>      dvv = _p[0][2];
      E = du * du;
      F = du * dv;
      G = dv * dv;
      e = N  * duu;
      f = N  * duv;
      g = N  * dvv;
  }


 } // END namespace GMlib
