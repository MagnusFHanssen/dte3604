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
}
