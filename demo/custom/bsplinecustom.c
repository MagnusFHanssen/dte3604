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

    int i;
    if (t == this->getEndP()){
        i = _t.size() - 4;
    } else{
        i = std::distance(_t.begin(), std::upper_bound(_t.begin(), _t.end(), t)) -1;
    }
    Vector<T,3> b = B(t, i);

    this->_p[0] = b[0] * _c[i-2]
                + b[1] * _c[i-1]
                + b[2] * _c[i];
    std::cout << this->_p[0] << std::endl;

}

template <typename T>
T BSplineCustom<T>::getStartP() const {return T(0);}

template <typename T>
T BSplineCustom<T>::getEndP() const {return T(1);}

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
void BSplineCustom<T>::makeKnots(int n){
    _t = {T(0), T(0), T(0)};
    T frac = T(1)/T(n-2);

    for (int i = 1; i < n-2; i++){
        _t.push_back(T(i*frac));
    }

    for (int i = 0; i < 3; i++){
        _t.push_back(T(1));
    }
}

}


