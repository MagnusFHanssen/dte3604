#include "blendingsurface.h"

namespace Custom{

using namespace GMlib;

template <typename T> inline
BlendingSurface<T>::BlendingSurface(PSurf<T,3>* copy, int u, int v)
    : PCurve<T,3>(20, 20, 1, 1), _startU{copy->getStartPU()},
    _endU{copy->getEndPU()}, _startV{copy->getStartPV()},
    _endV{copy->getEndPV()}, _closedU{copy->isClosedU()},
    _closedV{copy->isClosedV()}
{
    _fillKnotVector(_tU, _startU, _endU, u, _closedU);
    _fillKnotVector(_tV, _startV, _endV, v, _closedV);


}

template <typename T>
BlendingSurface<T>::~BlendingSurface(){}

// Getters
template <typename T>
T BlendingSurface<T>::getStartPU() const {return _startU;}

template <typename T>
T BlendingSurface<T>::getEndPU() const {return _endU;}

template <typename T>
T BlendingSurface<T>::getStartPV() const {return _startV;}

template <typename T>
T BlendingSurface<T>::getEndPV() const {return _endV;}

template <typename T>
bool BlendingSurface<T>::isClosedU() const {return _closedU;}

template <typename T>
bool BlendingSurface<T>::isClosedV() const {return _closedV;}

template <typename T>
void BlendingSurface<T>::_fillKnotVector(std::vector<T>& t, const T startP, const T endP, int steps, bool closed){
    t = {startP, startP};

    int n = steps - (closed ? 0 : 1);

    T frac = (endP - startP) / n;

    for (int i = 1; i < n; i++){
        t.push_back(startP + i * frac);
    }

    t.push_back(endP);
    t.push_back(endP);

    if (closed){
        t[0] = startP - frac;
        t[n + 2] = end + frac;
    }
}

}
