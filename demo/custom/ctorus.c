#include "ctorus.h"

namespace Custom{

using namespace GMlib;

template <typename T> inline
    CTorus<T>::CTorus(T a, T b) : PSurf<T, 3>(), _a{a}, _b{b}{}

template <typename T>
CTorus<T>::~CTorus(){}

template <typename T>
bool CTorus<T>::isClosedU() const {return true;}

template <typename T>
bool CTorus<T>::isClosedV() const {return true;}

template <typename T>
T CTorus<T>::getStartPU() const {return 0;}

template <typename T>
T CTorus<T>::getEndPU() const {return M_2PI;}

template <typename T>
T CTorus<T>::getStartPV() const {return 0;}

template <typename T>
T CTorus<T>::getEndPV() const {return M_2PI;}

template <typename T>
void CTorus<T>::eval( T u, T v, int d1, int d2, bool lu, bool lv) const{
    this->_p.setDim(d1 + 1, d2 + 2);

    this->_p[0][0][0] = (_b + _a * cos(v)) * cos(u);
    this->_p[0][0][1] = (_b + _a * cos(v)) * sin(u);
    this->_p[0][0][2] = _a * sin(v);

    if(this->_dm == GM_DERIVATION_EXPLICIT){
        if(d1){ // Su
            this->_p[1][0][0] =  -(_b + _a * cos(v)) * sin(u);
            this->_p[1][0][1] =  (_b + _a * cos(v)) * cos(u);
            this->_p[1][0][2] =  0;
        }
        if(d2){ // Sv
            this->_p[0][1][0] =  -_a * sin(v) * cos(u);
            this->_p[0][1][1] =  -_a * sin(v) * sin(u);
            this->_p[0][1][2] =  _a * cos(v);
        }
    }
}
}
