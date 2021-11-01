/**********************************************************************************
**
** Copyright (C) 1994 - 2016 University of Tromsø - The Arctic University of Norway
** Contact: GMlib Online Portal at https://source.uit.no/gmlib/gmlib/wikis/home
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





namespace GMlib {


//*****************************************
// Constructors and destructor           **
//*****************************************

/*! PCircle<T>::PCircle( T radius )
 *  Default constructor, to make a circle with the radius = radius.
 *
 *  \param[in] radius      (default 20) The radius of the circle
 */
  template <typename T>
  inline
  PCircle<T>::PCircle( T radius ) : PCurve<T,3>(20, 0, 7) {
      // Note that the last parameter in the PCurve constructor is 7,
      // this because 7 derivatives in eval() is implemented!!!!
    _rx = radius;
    _ry = radius;
  }



  /*! PCircle<T>::PCircle(  T radius1, T radius2 )
   *  Constructor, to make a ellipse with the horisontal radius = radius1
   *  and vertical radius = radius2
   *
   *  \param[in] radius      (default 20) The radius of the circle
   */
    template <typename T>
    inline
    PCircle<T>::PCircle(  T radius1, T radius2 ) : PCurve<T,3>(20, 0, 7) {
        // Note that the last parameter in the PCurve constructor is 7,
        // this because 7 derivatives in eval() is implemented!!!!
      _rx = radius1;
      _ry = radius2;
    }



  /*! PCircle<T>::PCircle(const PCircle<T>& copy )
   *  A copy constructor
   *  Making a copy of the curve (circle)
   *
   *  \param[in] copy The curve to copy
   */
  template <typename T>
  inline
  PCircle<T>::PCircle( const PCircle<T>& copy ) : PCurve<T,3>(copy) {}



  /*! PCircle<T>::~PCircle()
   *  The destructor
   *  clean up and destroy all private data
   */
  template <typename T>
  PCircle<T>::~PCircle() {}


  //**************************************
  //        Public local functons       **
  //**************************************

  /*! T PCircle<T>::getRadius() const
   *  Give you the radius of the circle (curve)
   *
   *  \return  The radius of the circle
   */
  template <typename T>
  inline
  T PCircle<T>::getRadius(int i) const {
    if(i==1)
      return _rx;
    else {
      return _ry;
    }
  }



  /*! void PCircle<T>::setRadius( T radius )
   *  Will change the radius of the circle
   *
   *  \param[in] radius   (default 20) The new radius of the circle
   */
  template <typename T>
  inline
  void PCircle<T>::setRadius( T radius ) {
      _rx = _ry = radius;
  }




    /*! void PCircle<T>::setRadius(T radius1, T radius2)
     *  Will change the radius of the ellipse
     *
     *  \param[in] radius   (default 20) The new radius of the circle
     */
    template <typename T>
    inline
    void PCircle<T>::setRadius( T radius1, T radius2 ) {
        _rx = radius1;
        _ry = radius2;
    }



  //***************************************************
  // Overrided (public) virtual functons from PCurve **
  //***************************************************

  /*! bool PCircle<T>::isClosed() const
   *  To tell that this curve (circle) is closed.
   */
  template <typename T>
  bool PCircle<T>::isClosed() const {
    return true;
  }



  //******************************************************
  // Overrided (protected) virtual functons from PCurve **
  //******************************************************

  /*! void PCircle<T>::eval( T t, int d, bool l ) const
   *  Evaluation of the curve at a given parameter value
   *  To compute position and d derivatives at parameter value t on the curve.
   *  7 derivatives are implemented
   *
   *  \param  t[in]  The parameter value to evaluate at
   *  \param  d[in]  The number of derivatives to compute
   *  \param  l[in]  (dummy) because left and right are always equal
   */
  template <typename T>
  void PCircle<T>::eval( T t, int d, bool /*l*/ ) const {

    this->_p.setDim( d + 1 );

    const T ct_x = _rx * cos(t);
    const T st_y = _ry * sin(t);

    this->_p[0][0] = ct_x;
    this->_p[0][1] = st_y;
    this->_p[0][2] = T(0);

    if( this->_dm == GM_DERIVATION_EXPLICIT ) {
      const T st_x = _rx * sin(t);
      const T ct_y = _ry * cos(t);
      if( d > 0 ) {
        this->_p[1][0] = -st_x;
        this->_p[1][1] =  ct_y;
        this->_p[1][2] =  T(0);
      }
      if( d > 1 ) this->_p[2] = -this->_p[0];
      if( d > 2 ) this->_p[3] = -this->_p[1];
      if( d > 3 ) this->_p[4] =  this->_p[0];
      if( d > 4 ) this->_p[5] =  this->_p[1];
      if( d > 5 ) this->_p[6] =  this->_p[2];
      if( d > 6 ) this->_p[7] =  this->_p[3];
    }
  }



  /*! T PCircle<T>::getStartP() const
   *  Provides the start parameter value associated with
   *  the eval() function implemented above.
   *  (the start parameter value = 0).
   *
   *  \return The parametervalue at start of the internal domain
   */
  template <typename T>
  T PCircle<T>::getStartP() const {
    return T(0);
  }



  /*! T PCircle<T>::getEndP() const
   *  Provides the end parameter value associated with
   *  the eval() function implemented above.
   *  (the end parameter value = 2*Pi).
   *
   *  \return The parametervalue at end of the internal domain
   */
  template <typename T>
  T PCircle<T>::getEndP()const {
    return T( M_2PI );
  }



  /*! void PCircle<T>::computeSurroundingSphere( const std::vector<DVector<Vector<T,3>>>&, Sphere<T,3>& s ) const
   *  Will set the surrounding sphere of the circle (equal the circle)
   *
   *  \param[in]  p   dummy!!!!, not used
   *  \param[out] s  The surrounding sphere to be computed
   */
  template <typename T>
  void PCircle<T>::computeSurroundingSphere( const std::vector<DVector<Vector<T,3>>>& /*p*/, Sphere<T,3>& s ) const {

     s.resetPos(Point<T,3>(T(0)));
     if(_rx > _ry)
       s.resetRadius(_rx);
     else
       s.resetRadius(_ry);
  }

} // END namespace GMlib
