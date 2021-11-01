#include "lotus.h"

#include <core/types/gmangle.h>

namespace GMlib{

Lotus::~Lotus(){}

Lotus::Lotus(double size):_s{size}{}

bool Lotus::isClosed() const {return true;}

double Lotus::getStartP() const {return 0.0;}
double Lotus::getEndP() const {return 1.0;}

void Lotus::eval(double t, int d, bool left) const
{
    this->_p.setDim(d + 1);

    const auto x = _s * (cos(16 * M_2PI * t) + cos(6 * M_2PI * t) / 2.0
                    + sin(10 * M_2PI * t) / 2.0);
    const auto y = _s * (sin(16 * M_2PI * t) + sin(6 * M_2PI * t) / 2.0
                    + cos(10 * M_2PI * t) / 2.0);

    this->_p[0][0] = x;
    this->_p[0][1] = y;
    this->_p[0][2] = 0.0;
}
}

