#ifndef MODELSURFACE_H
#define MODELSURFACE_H

#include <parametrics/gmpsurf.h>

namespace Custom {

using namespace GMlib;

template <typename T> class ModelSurface : public PSurf<T, 3>
{
    GM_SCENEOBJECT(ModelSurface)
public:
    ModelSurface();
    ~ModelSurface();

    bool isClosedU() const override;
    bool isClosedV() const override;

public:
    void eval( T u, T v, int d1, int d2, bool lu = true, bool lv = true ) const override;

    T getStartPU() const override;
    T getEndPU() const override;
    T getStartPV() const override;
    T getEndPV() const override;

    void init();

private:
    T _a;
    T _b;
    unsigned int _n;
};
}

#include "modelsurface.c"
#endif // MODELSURFACE_H
