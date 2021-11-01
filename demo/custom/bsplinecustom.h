#ifndef BSPLINECUSTOM_H
#define BSPLINECUSTOM_H

#include "../../gmlib/modules/parametrics/gmpcurve.h"

namespace GMlib {

    template<typename T>
    class BSplineCustom : public PCurve<T,3>
    {
        GM_SCENEOBJECT(BSplineCustom)
    public:
        BSplineCustom(const DVector<Point<T,3>>& c);
        BSplineCustom(const DVector<Point<T,3>>& p, int n);
        ~BSplineCustom();

        bool isClosed() const override;

    protected:
        void eval(T t, int d, bool left) const override;
        T getStartP() const override;
        T getEndP() const override;
    };
}
#endif // BSPLINECUSTOM_H
