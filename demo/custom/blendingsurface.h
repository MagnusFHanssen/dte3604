#ifndef BLENDINGSURFACE_H
#define BLENDINGSURFACE_H

#include <parametrics/gmpsurf.h>
#include "subpatch.h"

namespace Custom{

using namespace GMlib;

template <typename T>
class BlendingSurface : public PSurf<T, 3>
{
    GM_SCENEOBJECT(BlendingSurface)
public:
    BlendingSurface(PSurf<T, 3>* copy, int u, int v);
    ~BlendingSurface();

    bool isClosedU() const override;
    bool isClosedV() const override;

    void showLocalSurfaces();

protected:
    void eval( T u, T v, int d1, int d2, bool lu = true, bool lv = true ) const override;

    T getStartPU() const override;
    T getEndPU() const override;
    T getStartPV() const override;
    T getEndPV() const override;

    void init();

private:
    bool _closedU;
    bool _closedV;

    T _startU;
    T _endU;
    T _startV;
    T _endV;

    PSurf<T,3>* _surface;

    std::vector<T> _tU;
    std::vector<T> _tV;

    std::vector<std::vector<SubPatch<T>*>> _s;

    void _fillKnotVector(std::vector<T> &t, const T startP, const T endP, int steps, bool closed);

    void _makeLocalSurfaces(PSurf<T,3>* surface);

    T _getIndex(T t, bool is_u) const;

    T w(std::vector<T> const& knots, T const t, int i, int d) const;

    Vector<T,2> _blend(T w) const;

};
}

#include "blendingsurface.c"
#endif // BLENDINGSURFACE_H
