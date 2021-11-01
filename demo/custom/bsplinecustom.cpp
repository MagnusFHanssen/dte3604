#include "bsplinecustom.h"

namespace GMlib{

template <typename T>
BSplineCustom<T>::~BSplineCustom(){}

template <typename T>
BSplineCustom<T>::BSplineCustom(const DVector<Point<T,3>>& c)
{

}

template <typename T>
BSplineCustom<T>::BSplineCustom(const DVector<Point<T,3>>& p, int n)
{

}

template <typename T>
bool BSplineCustom<T>::isClosed() const {return false;}

template <typename T>
void BSplineCustom<T>::eval(T t, int d, bool left) const
{

}

template <typename T>
T BSplineCustom<T>::getStartP() const {return T(0);}

template <typename T>
T BSplineCustom<T>::getEndP() const {return T(1);}


}
