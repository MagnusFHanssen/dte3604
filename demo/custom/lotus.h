#ifndef LOTUS_H
#define LOTUS_H

//#include "../../gmlib/modules/parametrics/gmpcurve.h"
#include <parametrics/gmpcurve.h>

namespace GMlib {


class Lotus : public PCurve<double,3>
{
    GM_SCENEOBJECT(Lotus)
public:
    Lotus(double size=1.0);
    ~Lotus();

    bool isClosed() const override;

protected:
    void eval(double t, int d, bool left) const override;
    double getStartP() const override;
    double getEndP() const override;

private:
    double _s;

};


}


#endif // LOTUS_H
