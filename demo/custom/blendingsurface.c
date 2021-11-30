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

}
