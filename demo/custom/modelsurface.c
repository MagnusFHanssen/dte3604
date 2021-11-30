#include "modelsurface.h"

namespace Custom {

using namespace GMlib;

template <typename T> inline
    ModelSurface<T>::ModelSurface(unsigned int n, T a, T b) : PSurf<T, 3>(),
    _n{n}, _a{a}, _b{b}{}

template <typename T>
ModelSurface<T>::~ModelSurface(){}

template <typename T>
bool ModelSurface<T>::isClosedU() const {
    return true;
}

template <typename T>
bool ModelSurface<T>::isClosedV() const {
    return true;
}

template <typename T>
void ModelSurface<T>::eval( T u, T v, int d1, int d2, bool lu, bool lv) const{
    this->_p.setDim(d1 + 1, d2 + 1);

    T gamma = _a + _b * sin(2 * _n * v);
    T gamma_dv = 2 * _b * _n * cos(2 * _n * v);

    T x = cos(u + v) * cos(gamma);
    T x_du = -sin(u + v) * cos(gamma);
    T x_dv = x_du - cos(u + v) * sin(gamma) * gamma_dv;

    T y = sin(u + v) * cos(gamma);
    T y_du = x;
    T y_dv = y_du - sin(u + v) * sin(gamma) * gamma_dv;

    T z = cos(u - v) * cos(gamma);
    T z_du = -sin(u - v) * cos(gamma);
    T z_dv = sin(u - v) * cos(gamma) - cos(u - v) * sin(gamma) * gamma_dv;

    T w = sin(u - v) * sin(gamma);
    T w_du = cos(u - v) * sin(gamma);
    T w_dv = cos(u - v) * sin(gamma) + sin(u - v) * sin(gamma) * gamma_dv;

    T r = acos(w) / (M_PI * sqrt(1 - pow(w, 2)));
    T r_dw = - (1 - (w * acos(w))/sqrt(1 - pow(w, 2)))/(M_PI - M_PI * pow(w, 2));
    T r_du = w_du * r_dw;
    T r_dv = w_dv * r_dw;

    this->_p[0][0][0] = r * x;
    this->_p[0][0][1] = r * y;
    this->_p[0][0][2] = r * z;

    if(this->_dm == GM_DERIVATION_EXPLICIT){
        if(d1){ // Su
            this->_p[1][0][0] =  r_du * x + r * x_du;
            this->_p[1][0][1] =  r_du * y + r * y_du;
            this->_p[1][0][2] =  r_du * z + r * z_du;
        }
        if(d2){ // Sv
            this->_p[0][1][0] =  r_dv * x + r * x_dv;
            this->_p[0][1][1] =  r_dv * y + r * y_dv;
            this->_p[0][1][2] =  r_dv * z + r * z_dv;
        }
    }
}

template <typename T>
T ModelSurface<T>::getStartPU() const { return 0;}

template <typename T>
T ModelSurface<T>::getEndPU() const { return T(M_2PI);}

template <typename T>
T ModelSurface<T>::getStartPV() const { return 0;}

template <typename T>
T ModelSurface<T>::getEndPV() const { return T(M_PI);}

}
