#ifndef LANERIESENFELDCLOSED_H
#define LANERIESENFELDCLOSED_H

#include <parametrics/gmpcurve.h>

namespace Custom {

using namespace GMlib;

template <typename T> class LaneRiesenfeldClosed : public PCurve<T, 3> {
  GM_SCENEOBJECT(LaneRiesenfeldClosed)
public:
  LaneRiesenfeldClosed(const DVector<Point<T, 3>> &p, int d);
  ~LaneRiesenfeldClosed();

  bool isClosed() const override { return true; }

  void resample(std::vector<DVector<Vector<T, 3>>> &p, Sphere<T, 3> &s,
                const std::vector<T> &t, int d) const override;

protected:
  void eval(T t, int d, bool left) const override {}
  T getStartP() const override { return 0; }
  T getEndP() const override { return 1; }

private:
  DVector<Point<T, 3>> _p;
  int _d;

  void _doublePoints(std::vector<DVector<Vector<T, 3>>> &p,
                     int &no_elements) const;
  void _smoothPoints(std::vector<DVector<Vector<T, 3>>> &p, int d,
                     int no_elements) const;
};
} // namespace Custom
// Implementation
#include "laneriesenfeldclosed.c"
#endif // LANERIESENFELDCLOSED_H
