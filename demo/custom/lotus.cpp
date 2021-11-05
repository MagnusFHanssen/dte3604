#include "lotus.h"

#include <core/types/gmangle.h>

namespace Custom{

using namespace GMlib;

Lotus::~Lotus(){}

Lotus::Lotus(double size):PCurve<double,3>(300, 0, 0), _s{size}{}
Lotus::Lotus(const Lotus& copy):PCurve<double,3>(300, 0, 0), _s{copy._s}{}

bool Lotus::isClosed() const {return true;}

double Lotus::getStartP() const {return 0.0;}
double Lotus::getEndP() const {return 1.0;}

void Lotus::eval(double t, int d, bool left) const
{
    this->_p.setDim(d + 1);

    const auto y = _s * (cos(8 * M_2PI * t) + cos(3 * M_2PI * t) / 2.0
                    + sin(5 * M_2PI * t) / 3.0);
    const auto x = _s * (sin(8 * M_2PI * t) + sin(3 * M_2PI * t) / 2.0
                    + cos(5 * M_2PI * t) / 3.0);

    const auto z = _s * 0.2 * sin(9 * M_2PI * t);

    this->_p[0][0] = x;
    this->_p[0][1] = y;
    this->_p[0][2] = z;
}
}

