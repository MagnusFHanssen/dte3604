#ifndef BLENDINGSPLINE_H
#define BLENDINGSPLINE_H

#include <parametrics/curves/gmpsubcurve.h>
#include <parametrics/gmpcurve.h>

namespace Custom {

using namespace GMlib;

template <typename T> class BlendingSpline : public PCurve<T, 3> {
  GM_SCENEOBJECT(BlendingSpline)
public:
  BlendingSpline(PCurve<T, 3>* copy, int n);
  ~BlendingSpline();

  bool isClosed() const override;

  void showControlCurves();

protected:
  void eval(T t, int d, bool left) const override;
  T getStartP() const override;
  T getEndP() const override;

private:
  std::vector<PSubCurve<T> *> _c;
  std::vector<T> _t;
  bool _closed;

  void _makeKnots(int n, T start = 0, T end = 1);

  void _makeControlCurves(PCurve<T, 3>* curve);

  T _t_start;
  T _t_end;

  int _getIndex(T t) const;

  T w(T t, int i, int d) const;
  Vector<T, 2> B(T t, int i) const;

  T _blend(T w) const;

  Point<T,3> c(int i, T t) const;
};
} // namespace Custom
#include "blendingspline.c"
#endif // BLENDINGSPLINE_H
