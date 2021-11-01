/**********************************************************************************
**
** Copyright (C) 1994 Narvik University College
** Contact: GMlib Online Portal at http://episteme.hin.
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

namespace GMlib {


  //*****************************************
  // Constructors and destructor           **
  //*****************************************

  template <typename T>
  inline
  PBSplineSurf<T>::PBSplineSurf( const DMatrix< Vector<T,3> >& c, const DVector<T>& u,  const DVector<T>& v, int du, int dv) {

      init();
      _c = c;                           // Set Control Points
      _ku = u.getDim()-_c.getDim1();    // Set order
      _kv = v.getDim()-_c.getDim2();

      // Set knot-vector in u-direction
      if(_ku == 1) {  // closed in u-direction
          _du = du;
          initKnot( _u, _cu, _ku, u, _c.getDim1(), _du);
      } else {
          _u = u;
          _du = _ku-1;
          _cu = false;
      }
      // Set knot-vector in v-direction
      if(_kv == 1) {  // closed in v-direction
          if(dv==0)   _dv = du;
          else        _dv = dv;
          initKnot( _v, _cv, _kv, v, _c.getDim2(), _dv);
      } else {
          _v = v;
          _dv = _kv-1;
          _cv = false;
      }
  }


  template <typename T>
  inline
  PBSplineSurf<T>::PBSplineSurf( const PBSplineSurf<T>& copy ) : PSurf<T,3>( copy ) {

      init();

      _c = copy._c;
      _u = copy._u;
      _v = copy._v;
      _ku = copy._ku;
      _kv = copy._kv;
      _du = copy._du;
      _dv = copy._dv;
      _cu = copy._cu;
      _cv = copy._cv;
  }


  template <typename T>
  PBSplineSurf<T>::~PBSplineSurf() {

      if(_sgv) delete _sgv;
  }


  //*********************************
  //**   Public local functons     **
  //*********************************

  template <typename T>
  inline
  DMatrix< Vector<T,3> >& PBSplineSurf<T>::getControlPoints() const {
      return _c;
  }


  template <typename T>
  inline
  int PBSplineSurf<T>::getDegreeU() const {
      return _du;
  }


  template <typename T>
  inline
  int PBSplineSurf<T>::getDegreeV() const {
      return _dv;
  }


  template <typename T>
  bool PBSplineSurf<T>::isSelectorsVisible() const {
      return _selectors;
  }



  template <typename T>
  void PBSplineSurf<T>::setClosed( bool closed_u, bool closed_v, T du, T dv ) {

      bool changed = false;

      if(_cu != closed_u) {
          DVector<T> nu = _u;
          if(closed_u) {
              if(du==T(0)) du = (_u[_u.getDim()-_ku]-_u[_du])/getIntervallInVector(_u);
              initKnot( _u, _cu, _ku, nu, _c.getDim1(), _du, du);
          } else
              initKnot2(_u, _cu, nu, _c.getDim1(), _du);
          changed = true;
      }

      if(_cv != closed_v){
          DVector<T> nv = _v;
          if(closed_v) {
              if(dv==T(0)) dv = (_v[_v.getDim()-_kv]-_v[_dv])/getIntervallInVector(_v);
              initKnot( _v, _cv, _kv, nv, _c.getDim2(), _dv, dv);
          } else
              initKnot2(_v, _cv, nv, _c.getDim2(), _dv);
          changed = true;
      }

      if(changed && _selectors) {
          hideSelectors();
          showSelectors(_selector_radius, _grid, _selector_color, _grid_color);
      }
  }



  template <typename T>
  inline
  void PBSplineSurf<T>::setPartitionCritere(int pcu, int pcv) {
      _pcu = pcu;
      _pcv = pcv;
  }



  template <typename T>
  inline
  void PBSplineSurf<T>::setControlPoints( const DMatrix< Vector<T,3> >& cp ) {

      if( _c.getDim1() == cp.getDim1() && _c.getDim2() == cp.getDim2() )
          _c = cp;
      else
          std::cerr << "Can not change the control point because the dimentions are wrong!";
  }



  template <typename T>
  inline
  void PBSplineSurf<T>::updateCoeffs( const Vector<T,3>& d ) {

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
          this->translateParent( -d );
      }
  }



  template <typename T>
  inline
  void PBSplineSurf<T>::enablePartitionVisualizer(int u, int v) {
      _part_viz = true;
      setPartitionCritere(u, v);
      this->enableDefaultVisualizer(true);
}



  //********************************************************
  // Overrided (public) virtual functons from SceneObject **
  //********************************************************

  template <typename T>
  void PBSplineSurf<T>::edit( int selector_id, const Vector<T,3>& dp  ) {

    _c_moved = true;
      if( this->_parent ) this->_parent->edit( this );
      _pos_change.push_back(EditSet(selector_id, dp));
      this->setEditDone();
    _c_moved = false;
  }



//  template <typename T>
//  inline
//  void PBSplineSurf<T>::insertVisualizer( Visualizer *visualizer ) {

//    SceneObject::insertVisualizer( visualizer );

//    PSurfVisualizer<T,3> *visu = dynamic_cast<PSurfVisualizer<T,3>*>( visualizer );
//    if( !visu )
//      return;

//    for(int i=0; i < _visuS.getDim1(); i++)
//        for(int j=0; j < _visuS.getDim2(); j++) {
//            _visuS[i][j].vis += visu;
//        }
//  }




  template <typename T>
  void PBSplineSurf<T>::replot( ) const{

    updateSamples();
    PSurf<T,3>::replot();
  }




  //**************************************************
  // Overrided (public) virtual functons from PSurf **
  //**************************************************

  template <typename T>
  void PBSplineSurf<T>::sample( int m1, int m2, int d1, int d2 ) {

    if(_part_viz) {
      this->cleanVisualizers();
      initSample(m1, m2, d1, d2);
      makePartition();
      for(int i=0; i < _visu.getDim1(); i++)
        for(int j=0; j < _visu.getDim2(); j++) {
          DMatrix<DMatrix<Vector<T,3>>>&  p = _visu[i][j].sample_val;
          DMatrix<Vector<float,3>>& normals = _visu[i][j].normals;
          Sphere<T,3>&                    s = _visu[i][j].sur_sphere;
          // Calculate sample positions, derivatives, normals and surrounding sphere
          resample(p, _pre_basis_u[i], _pre_basis_v[j], _vpu[i].m, _vpv[j].m, d1, d2 );
          this->resampleNormals( p, normals );
          this->computeSurroundingSphere(p, s);
        }
    } else
      PSurf<T,3>::sample( m1, m2, d1, d2 );
    preSample();
  }




  /*! void  PBSplineSurf<T>::preSample()
   *  Private, not for public use
   *  Compute all sample points for all partitions
   *
   *  \param[in]  d   The number of derivatives to compute
   */
  template <typename T>
  void  PBSplineSurf<T>::preSample() {

    _cp_index.setDim(_c.getDim1(), _c.getDim2());
    for(int i=0; i < _cp_index.getDim1(); i++)
      for(int j=0; j < _cp_index.getDim2(); j++) {
        _cp_index[i][j].setDim(_pre_basis_u.size(), _pre_basis_v.size());
        for(int i1=0; i1 < _cp_index[i][j].getDim1(); i1++)
          for(int j1=0; j1 < _cp_index[i][j].getDim2(); j1++) {
            _cp_index[i][j][i1][j1][0] = {std::numeric_limits<int>::max(),-1};
            _cp_index[i][j][i1][j1][1] = {std::numeric_limits<int>::max(),-1};
          }
      }
    for(uint i=0; i<_pre_basis_u.size(); i++)
      for(uint j=0; j<_pre_basis_v.size(); j++)
        for(uint i1=0; i1<_pre_basis_u[i].size(); i1++)
          for(uint j1=0; j1<_pre_basis_v[j].size(); j1++)
            for(uint k1=0; k1<_pre_basis_u[i][i1].ind.size(); k1++ )
              for(uint k2=0; k2<_pre_basis_v[j][j1].ind.size(); k2++ ) {
                int i_p = _pre_basis_u[i][i1].ind[k1];
                int j_p = _pre_basis_v[j][j1].ind[k2];
                if(int(i1) < _cp_index[i_p][j_p][i][j][0][0]) _cp_index[i_p][j_p][i][j][0][0] = i1;
                if(int(i1) > _cp_index[i_p][j_p][i][j][0][1]) _cp_index[i_p][j_p][i][j][0][1] = i1;
                if(int(j1) < _cp_index[i_p][j_p][i][j][1][0]) _cp_index[i_p][j_p][i][j][1][0] = j1;
                if(int(j1) > _cp_index[i_p][j_p][i][j][1][1]) _cp_index[i_p][j_p][i][j][1][1] = j1;
              }
  }



  template <typename T>
  bool PBSplineSurf<T>::isClosedU() const {
      return _cu;
  }


  template <typename T>
  bool PBSplineSurf<T>::isClosedV() const {
      return _cv;
  }


  template <typename T>
  void PBSplineSurf<T>::showSelectors( T rad, bool grid, const Color& selector_color, const Color& grid_color ) {

      if(! _selectors ) {
          _s.setDim( _c.getDim1(), _c.getDim2() );

          for( int i = 0; i < _c.getDim1(); i++ )
              for( int j = 0; j < _c.getDim2(); j++ )
                  this->insert( _s[i][j] = new Selector<T,3>( _c[i][j], _map2(i, j), this, rad, selector_color ) );
          _selectors       = true;
      }
      _selector_radius = rad;
      _selector_color  = selector_color;

      if( grid ) {
          if(!_sgv) _sgv = new SelectorGridVisualizer<T>;
          _sgv->setSelectors( _c, _cu, _cv );
          _sgv->setColor( grid_color );
          SceneObject::insertVisualizer( _sgv );
          this->setEditDone();
      }
      _grid_color      = grid_color;
      _grid            = grid;
  }



  template <typename T>
  void PBSplineSurf<T>::hideSelectors() {

      // Remove Selectors
      if( _selectors ) {
          for( int i = 0; i < _s.getDim1(); i++ ) {
              for( int j = 0; j < _s.getDim2(); j++ ) {
                  this->remove( _s[i][j] );
                  delete _s[i][j];
              }
          }
          _s.resetDim(0,0);
          _selectors = false;
      }

      // Hide Selector Grid Visualizer
      if(_sgv){
        this->removeVisualizer( _sgv );
        _sgv->reset();
      }
  }




  /*! void PBSplineSurf<T>::toggleSelectors()
   *  To toggle all selectors and selector grid.
   */
  template <typename T>
  void PBSplineSurf<T>::toggleSelectors() {

    if(_selectors)  hideSelectors();
    else            showSelectors();
  }





  /*! void PBSplineSurf<T>::toggleClose()
   *  To toggle open surface - closed surface
   */
  template <typename T>
  void PBSplineSurf<T>::toggleClose() {

      if(!isClosedU() && _ku >_c.getDim1()) return;
      if(!isClosedV() && _kv >_c.getDim2()) return;

      setClosed(!isClosedU(), !isClosedV() );
      preSample(1, _visu.no_sample[0]);
      preSample(2, _visu.no_sample[1]);
      sample( _visu.no_sample[0], _visu.no_sample[1], _visu.no_derivatives[0], _visu.no_derivatives[1] );
      if(_selectors) {
          hideSelectors();
          showSelectors(_selector_radius, _grid, _selector_color, _grid_color);
      }
      this->setEditDone();
  }



  //*****************************************************
  // Overrided (protected) virtual functons from PSurf **
  //*****************************************************


  template <typename T>
  void PBSplineSurf<T>::eval( T u, T v, int du, int dv, bool lu, bool lv ) const {

      DMatrix<T>   bu, bv;
      std::vector<int> ind_i(_ku), ind_j(_kv);

      int i = EvaluatorStatic<T>::evaluateBSp( bu, u, _u, _du, lu) - _du;
      int j = EvaluatorStatic<T>::evaluateBSp( bv, v, _v, _dv, lv) - _dv;

      makeIndex(ind_i, i, _ku, _c.getDim1());
      makeIndex(ind_j, j, _kv, _c.getDim2());

      multEval( this->_p, bu, bv, ind_i, ind_j, du, dv);
  }



  template <typename T>
  T PBSplineSurf<T>::getEndPU() const {
      return _u(_u.getDim()-_ku);
  }



  template <typename T>
  T PBSplineSurf<T>::getEndPV() const {
      return _v(_v.getDim()-_kv);
  }



  template <typename T>
  T PBSplineSurf<T>::getStartPU() const {
      return _u(_du);
  }



  template <typename T>
  inline
  T PBSplineSurf<T>::getStartPV() const {
      return _v(_dv);
  }



  //*****************************************
  //     Local (protected) functons        **
  //*****************************************


  template <typename T>
  inline
  void PBSplineSurf<T>::init() {

      this->_type_id = GM_SO_TYPE_SURFACE_BSPLINE;

      _selectors    = false;
      _c_moved      = false;

      _part_viz       = false;
      _pcu = 1;
      _pcv = 1;

      _pre_basis_u.resize(1);
      _pre_basis_u.resize(1);

      _sgv =  0x0;
  }




  //***************************************************
  // Overrided (private) virtual functons from PSurf **
  //***************************************************


  // Sampling for pre-evaluation of basis functions, chose direction
  //*****************************************************************
  template <typename T>
  void  PBSplineSurf<T>::preSample( int dir, int m ) {

    if(_part_viz) {
      if( dir==1 ) {
        preparePartition( _vpu, _u, _ku, _pcu, m );
        _pre_basis_u.resize(_vpu.getDim());
        for(int i=0; i<_pre_basis_u.size(); i++)
          preSample( _pre_basis_u[i], _u, _vpu[i].m, _du, _c.getDim1(), _u[_vpu[i].is], _u[_vpu[i].ie] );
      }
      if( dir==2 ) {
        preparePartition( _vpv, _v, _kv, _pcv, m );
        _pre_basis_v.resize(_vpv.getDim());
        for(int i=0; i<_pre_basis_v.size(); i++)
          preSample( _pre_basis_v[i], _v, _vpv[i].m, _dv, _c.getDim2(), _v[_vpv[i].is], _v[_vpv[i].ie] );
      }
    }
    else {
      if( dir==1 ) {
        _pre_basis_u.resize(1);
        preSample( _pre_basis_u[0], _u, m, _du, _c.getDim1(), _u[_du], _u[_u.getDim()-_ku] );
      }
      if( dir==2 ) {
        _pre_basis_v.resize(1);
        preSample( _pre_basis_v[0], _v, m, _dv, _c.getDim2(), _v[_dv], _v[_v.getDim()-_kv] );
      }
    }
  }





  template <typename T>
  void PBSplineSurf<T>::resample( DMatrix< DMatrix< Vector<T,3> > >& p, const PreBasis<T>& bu, const PreBasis<T>& bv, int m1, int m2, int d1, int d2 ) const {

      p.setDim(m1, m2);

      for(int i=0; i<m1; i++)
          for(int j=0; j<m2; j++)
              multEval( p[i][j], bu[i], bv[j], bu[i].ind, bv[j].ind, d1, d2);
  }




  //***************************************
  //**   Local (private) help functons   **
  //***************************************

  template <typename T>
  inline
  void PBSplineSurf<T>::makeIndex( std::vector<int>& ind, int i, int k, int n) const{

      if(i+k > n){
          int j, s = n-i;
          for (j=0; j < s; j++)
              ind[j] = i++;
          for ( ; j < k; j++)
              ind[j] = j-s;
      }
      else
          for(int j=0; j<k; j++)
              ind[j]= i+j;
  }



  template <typename T>
  inline
  void PBSplineSurf<T>::multEval(DMatrix<Vector<T,3>>& p, const DMatrix<T>& bu, const DMatrix<T>& bv, const std::vector<int>& ii, const std::vector<int>& ij, int du, int dv) const {

      // Set Dimensions
      p.setDim(du+1,dv+1);

      // We manually do these two operations!
      //    bv.transpose();
      //    p = bu * (c^bv);
      // The reason why we do it manually is that we only copmpute the du first lines of bu and the dv first lines of bv

      //    c = _c^bvT
      DMatrix<Vector<T,3>> c(_ku, dv+1);

      for(int i=0; i< _ku; i++)
          for(int j=0; j<=dv; j++) {
              c[i][j] = _c(ii[i])(ij[0])*bv(j)(0);
              for(int k=1; k<_kv; k++)
                  c[i][j] += _c(ii[i])(ij[k])*bv(j)(k);
          }
      //    p = bu * c
      for(int i=0; i<=dv; i++)
          for(int j=0; j<=du; j++) {
              p[i][j] = bu(i)(0)*c[0][j];
              for(int k=1; k<_ku; k++)
                  p[i][j] += bu(i)(k)*c[k][j];
          }
  }




  template <typename T>
  inline
  void PBSplineSurf<T>::initKnot( DVector<T>& t, bool& c, int& k, const DVector<T>& g, int n, int d, T dt) {

      k = d+1;
      t.setDim(n+d+k);
      if(dt > T(0)) {
          for(int i = d; i < n+1; i++)                  t[i] = g(i);
          for(int i = n+1; i < n+k; i++)                t[i] = t[i-1] + dt;
      } else
          for(int i = d; i < n+k; i++)                  t[i] = g(i-d);
      for(int i = d-1, j = n+d; i >= 0; i--,j--)        t[i] = t[i+1] - (t[j]-t[j-1]);
      for(int i = n+k, j = k; i < t.getDim(); i++,j++)  t[i] = t[i-1] + (t[j]-t[j-1]);
      c = true;
  }




  template <typename T>
  inline
  void PBSplineSurf<T>::initKnot2( DVector<T>& t, bool& c, const DVector<T>& g, int n, int d) {

      t.setDim(g.getDim()-d);
      for(int i = d; i <= n; i++)            t[i] = g(i);
      for(int i = d-1; i >= 0; i--)          t[i] = t[i+1];
      for(int i = n+1; i < t.getDim(); i++)  t[i] = t[i-1];
      c = false;
  }




  // pre-evaluation of basis fuctions, the B-spline-Hermite matrix, independent of direction
  //****************************************************************************************
  template <typename T>
  inline
  void PBSplineSurf<T>::preSample( PreBasis<T>& p, const DVector<T>& t, int m, int d, int n, T start, T end ) {

    const T dt = ( end - start ) / T(m-1); // dt is the step in parameter values
      p.resize(m);                         // p is a vector of  Bernstein-Hermite matrises at sample points

      // Compute the Bernstein-Hermite matrix
      for( int j = 0; j < m-1; j++ ) {
          int i = EvaluatorStatic<T>::evaluateBSp( p[j], start+j*dt, t, d, false );// - d;
          p[j].ind.init( i, d+1, n);
      }
      int i = EvaluatorStatic<T>::evaluateBSp( p[m-1], end, t, d, true );// - d;
      p[m-1].ind.init( i, d+1, n);
  }





  /*! void  PBSplineSurf<T>::updateSamples() const
   *  Private, not for public use
   *  Update affected sample points for all partitions when control points have been moved
   */
  template <typename T>
  void  PBSplineSurf<T>::updateSamples() const {

    while(_pos_change.size()>0) {
      EditSet es = _pos_change.back();
      for(int i=_pos_change.size()-2; i>=0; i--)
        if (_pos_change[i].ind == es.ind) {
          es.dp += _pos_change[i].dp;
          _pos_change.erase(_pos_change.begin()+i);
        }
      Vector<int,2> ind = _map1(es.ind);
      for(int i=0; i<_visu.getDim1(); i++)
        for(int j=0; j<_visu.getDim2(); j++) {
          // for each partition
          Partition& p = this->_visu[i][j];
          for(int i1 =_cp_index[ind[0]][ind[1]][i][j][0][0]; i1 <= _cp_index[ind[0]][ind[1]][i][j][0][1]; i1++)
            for(int j1 =_cp_index[ind[0]][ind[1]][i][j][1][0]; j1 <= _cp_index[ind[0]][ind[1]][i][j][1][1]; j1++) {
              // for each sample value
              DMatrix<Vector<T,3>>&  sp = p.sample_val[i1][j1];
              for(uint ku=0; ku < _pre_basis_u[i][i1].ind.size(); ku++)
                for(uint kv=0; kv < _pre_basis_v[j][j1].ind.size(); kv++)
                  if(_pre_basis_u[i][i1].ind[ku] == ind[0] && _pre_basis_v[j][j1].ind[kv] == ind[1])
                    comp(sp, _pre_basis_u[i][i1],  _pre_basis_v[j][j1], es.dp, ku, kv);
              p.normals[i1][j1] = sp(1)(0) ^ sp(0)(1);
              p.normals[i1][j1].normalize();
              p.sur_sphere += sp(0)(0);
            }
        }
      _pos_change.pop_back();
    }
  }



  /*! void PBSplineSurf<T>::comp(DVector<Vector<T,3>>& p, const DMatrix<T>& B, const Vector<T,3>& d, int k) const
   *  Protected,
   *  Partial vector-matrix computation, ie. actually a vector-vector innerproduct.
   *  where vi compute a vector vith a column vector with index ku in the matrix Bu,
   *  i.e. p = Bu^T_row[ku] * d * Bv_column[kv]
   *
   *  \param[out]  p  Updating the position and d derivatives
   *  \param[in]   Bu The Bernstein-Hermite matrix in u - direction.
   *  \param[in]   Bv The Bernstein-Hermite matrix in v - direction.
   *  \param[in]   d  The distance vector, one control point has been moved
   *  \param[in]   ku the u-index of the column to compute with, also a mapped index of the control point that has been moved
   *  \param[in]   ku the v-index of the column to compute with, also a mapped index of the control point that has been moved
   */
  template <typename T>
  inline
  void  PBSplineSurf<T>::comp(DMatrix<Vector<T,3>>& p, const DMatrix<T>& Bu, const DMatrix<T>& Bv, const Vector<T,3>& d, int ku, int kv) const {

      for(int i=0; i<p.getDim1(); i++)
        for(int j=0; j<p.getDim2(); j++)
          p[i][j] += (Bu(i)(ku)*Bv(j)(kv))*d;
  }





  template <typename T>
  inline
  void  PBSplineSurf<T>::makePartition() {

    // Inserting vizualizers
    _visu.setDim(_vpu.getDim(), _vpv.getDim());
    for(int i=0; i<_visu.getDim1(); i++)
      for(int j=0; j<_visu.getDim2(); j++) {
        _visu[i][j] = {_vpu[i].m,_vpv[j].m};
        _visu[i][j].s_e_u = {_u[_vpu[i].is], _u[_vpu[i].ie]};
        _visu[i][j].s_e_v = {_v[_vpv[j].is], _v[_vpv[j].ie]};
      }
    this->prepareVisualizers();
  }




  template <typename T>
  inline
  Vector<int,2> PBSplineSurf<T>::_map1(int i) const {

    Point<int,2> r;
    r[0] = i/_c.getDim2();
    r[1] = i - r[0]*_c.getDim2();
    return r;
  }



  template <typename T>
  inline
  int PBSplineSurf<T>::_map2(int i, int j) const {

      return i*_c.getDim2()+j;
  }



} // END namespace GMlib


