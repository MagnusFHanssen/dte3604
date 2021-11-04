#include "laneriesenfeldclosed.h"

namespace Custom {

using namespace GMlib;

template <typename T>
LaneRiesenfeldClosed<T>::LaneRiesenfeldClosed(const DVector<Point<T, 3>> &p,
                                              int d)
    : PCurve<T, 3>(100, 0, 0), _p{p}, _d{d} {}

template <typename T> LaneRiesenfeldClosed<T>::~LaneRiesenfeldClosed() {}

template <typename T>
void LaneRiesenfeldClosed<T>::resample(std::vector<DVector<Vector<T, 3>>> &p,
                                       Sphere<T, 3> &s, const std::vector<T> &t,
                                       int d) const {
  int n = _p.getDim();
  int k = this->_visu.no_sample;
  int m = n * (2 << k - 1) + 1;

  p.resize(m);
  s.reset();

  for (int i = 0; i < m; i++){
      p[i].setDim(1);
  }

  for (int i = 0; i < n; i++) {
    p[i][0] = _p[i];
  }
  p[n][0] = _p[0];

  for (int i = 0; i < k; i++) {
    _doublePoints(p, n);
    _smoothPoints(p, _d, n);
  }

  computeSurroundingSphere(p, s);
  if (d > _der_implemented || (d > 0 && this->_dm == GM_DERIVATION_DD))
    DD::compute1D(p, t, isClosed(), d, _der_implemented);
}

template <typename T>
void LaneRiesenfeldClosed<T>::_doublePoints(
    std::vector<DVector<Vector<T, 3>>> &p, int &no_elements) const {
  for (int i = no_elements; i > 0; --i) {
    p[2 * i][0] = p[i][0];
    p[2 * i - 1][0] = T(0.5) * (p[i][0] + p[i - 1][0]);
  }
  no_elements = 2 * no_elements;
}

template <typename T>
void LaneRiesenfeldClosed<T>::_smoothPoints(
    std::vector<DVector<Vector<T, 3>>> &p, int d, int no_elements) const {
  for (int j = 1; j < d; j++) {
    for (int i = 0; i < no_elements; i++) {
      p[i][0] = T(0.5) * (p[i][0] + p[i + 1][0]);
    }
    p[no_elements][0] = p[0][0];
  }
}

} // namespace Custom
