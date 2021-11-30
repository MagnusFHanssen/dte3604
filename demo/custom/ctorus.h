#ifndef CTORUS_H
#define CTORUS_H

#include <parametrics/gmpsurf.h>

namespace Custom{

using namespace GMlib;

template <typename T>
class CTorus : public PSurf<T, 3>
{
    GM_SCENEOBJECT(CTorus)
public:
    CTorus(T a = T(2), T b = T(3));
    ~CTorus();

    bool isClosedU() const override;
    bool isClosedV() const override;

public:
    void eval( T u, T v, int d1, int d2, bool lu = true, bool lv = true ) const override;

    T getStartPU() const override;
    T getEndPU() const override;
    T getStartPV() const override;
    T getEndPV() const override;

private:
    T _a;
    T _b;
};
}

#include "ctorus.c"
#endif // CTORUS_H
