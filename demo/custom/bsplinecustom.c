#include "bsplinecustom.h"

namespace Custom{

using namespace GMlib;

template <typename T>
BSplineCustom<T>::BSplineCustom(const DVector<Point<T,3>>& c, bool closed)
    :PCurve<T,3>(100, 0, 0),_c{c},_closed{closed},_useBlending{false}
{
    makeKnots(_c.getDim());
}

template <typename T>
BSplineCustom<T>::BSplineCustom(const DVector<Point<T,3>>& p, int n, bool closed)
    :PCurve<T,3>(100, 0, 0),_closed{closed},_useBlending{false}
{
    makeControlPoints(p, n);
}


template <typename T>
BSplineCustom<T>::BSplineCustom(const BSplineCustom<T>& copy):PCurve<T,3>(100, 0, 0),
    _c{copy._c},_t{copy._t},_closed{copy._closed},_useBlending{copy._useBlending}{}

template <typename T>
BSplineCustom<T>::~BSplineCustom(){}

template <typename T>
bool BSplineCustom<T>::isClosed() const {return _closed;}

template <typename T>
void BSplineCustom<T>::setBlending(bool b){_useBlending = b;}

template <typename T>
void BSplineCustom<T>::eval(T t, int d, bool left) const
{
    this->_p.setDim(d+1);

    int i = _getIndex(t);

    Vector<T,3> b = B(t, i);

    this->_p[0] = b[0] * _c[i-2]
                + b[1] * _c[i-1]
                + b[2] * _c[i];
}

template <typename T>
T BSplineCustom<T>::getStartP() const {return _t[2];}

template <typename T>
T BSplineCustom<T>::getEndP() const {return _t[_t.size()-3];}

template <typename T>
T BSplineCustom<T>::w(T t, int i, int d) const{
    return (t - _t[i])/(_t[i+d]-_t[i]);
}

template <typename T>
Vector<T,3> BSplineCustom<T>::B(T t, int i) const{
    T w_a, w_b, w_c;
    if (_useBlending){
        w_a = _blend(w(t, i, 1));
        w_b = _blend(w(t, i-1, 2));
        w_c = _blend(w(t, i, 2));
    }else{
        w_a = w(t, i, 1);
        w_b = w(t, i-1, 2);
        w_c = w(t, i, 2);
    }

    T b_1 = T(1) - w_a - w_b + w_a * w_b;
    T b_2 = w_b - w_a * w_b + w_a - w_a * w_c;
    T b_3 = w_a * w_c;

    return Vector<T,3>(b_1, b_2, b_3);
}

template <typename T>
T BSplineCustom<T>::_blend(T w) const {
    return w - T(1.0/(M_2PI))*sin(M_2PI*w);
}

template <typename T>
int BSplineCustom<T>::_getIndex(T t) const{
    int i;
    if (t == this->getEndP()){
        i = _t.size() - 4;
    } else{
        i = std::distance(_t.begin(), std::upper_bound(_t.begin(), _t.end(), t)) -1;
    }
    return i;
}



template <typename T>
void BSplineCustom<T>::makeKnots(int n, T start, T end){
    _t = {start, start, start};
    T frac = (end-start)/T(n-2);

    for (int i = 1; i < n-2; i++){
        _t.push_back(T(start + i*frac));
    }

    for (int i = 0; i < 3; i++){
        _t.push_back(end);
    }
}

template <typename T>
void BSplineCustom<T>::makeControlPoints(const DVector<Point<T,3>>& p, int n){
    std::vector<T> x = {0};
    for (int i = 1; i < p.getDim(); i++){
        x.push_back((p[i] - p[i-1]).getLength() + x.back());
    }

    makeKnots(n, x[0], x.back());

    DMatrix<T> A(p.getDim(), n, T(0));
    for (int i = 0; i < x.size(); i++){
        int j = _getIndex(x[i]);
        auto b = B(x[i], j);
        A[i][j-2] = b[0];
        A[i][j-1] = b[1];
        A[i][j]   = b[2];
    }
    auto A_t = A;

    A_t.transpose();

    auto B_mat = A_t * A;

    auto y = A_t * p;

    B_mat.invert();

    _c = B_mat * y;
}

}


