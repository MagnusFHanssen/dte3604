#ifndef BLENDINGSURFACE_H
#define BLENDINGSURFACE_H

#include <parametrics/gmpsurf.h>
#include <parametrics/surfaces/gmpsubsurf.h>

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




};
}
#endif // BLENDINGSURFACE_H
