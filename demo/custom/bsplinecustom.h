#ifndef BSPLINECUSTOM_H
#define BSPLINECUSTOM_H

#include "../../gmlib/modules/parametrics/gmpcurve.h"

namespace Custom {

using namespace GMlib;

template <typename T> class BSplineCustom : public PCurve<T, 3> {
  GM_SCENEOBJECT(BSplineCustom)
public:
  BSplineCustom(const DVector<Point<T, 3>> &c, bool closed = false);
  BSplineCustom(const DVector<Point<T, 3>> &p, int n, bool closed = false);
  BSplineCustom(const BSplineCustom<T> &copy);
  ~BSplineCustom();

  bool isClosed() const override;

  void setBlending(bool b = false);

protected:
  void eval(T t, int d, bool left) const override;
  T getStartP() const override;
  T getEndP() const override;

private:
  // Control points
  DVector<Point<T, 3>> _c;
  // Knot vector
  std::vector<T> _t;

  bool _closed;
  bool _useBlending;

  // Methods
  T w(T t, int i, int d) const;
  Vector<T, 3> B(T t, int i) const;

  T _blend(T w) const;

  int _getIndex(T t) const;

  void makeKnots(int n, T start = 0, T end = 1);
  void makeControlPoints(const DVector<Point<T, 3>> &p, int n);
};
} // namespace Custom
// Implementation
#include "bsplinecustom.c"

#endif // BSPLINECUSTOM_H
