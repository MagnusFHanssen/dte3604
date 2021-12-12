#include "blendingspline.h"

namespace Custom {

using namespace GMlib;

template <typename T>
BlendingSpline<T>::BlendingSpline(PCurve<T, 3> *copy, int n)
    : PCurve<T, 3>(20, 0, 0), _t_start{copy->getParStart()},
      _t_end{copy->getParEnd()}, _closed{copy->isClosed()} {

  _makeKnots(n, _t_start, _t_end);

  _makeControlCurves(copy);
}

template <typename T> BlendingSpline<T>::~BlendingSpline() {}

template <typename T> bool BlendingSpline<T>::isClosed() const {
  return _closed;
}

template <typename T> T BlendingSpline<T>::getStartP() const {
  return _t_start;
}

template <typename T> T BlendingSpline<T>::getEndP() const { return _t_end; }

template <typename T>
void BlendingSpline<T>::_makeKnots(int n, T start, T end) {
  _t = {start, start};

  int n_ = n - (_closed ? 0 : 1);

  T frac = (end - start) / n_;

  for (int i = 1; i < n_; i++) {
    _t.push_back(start + i * frac);
  }

  _t.push_back(end);
  _t.push_back(end);

  if (_closed) {
    _t[0] = start - frac;
    _t[n_ + 2] = end + frac;
  }
}

template <typename T>
void BlendingSpline<T>::_makeControlCurves(PCurve<T, 3> *curve) {
  int n = _t.size() - 2;
  _c.resize(n);
  if (_closed) { // If the curve is closed, the ends are the same
      n--;
    for (int i = 0; i < n; i++) {
      _c[i] = new PSubCurve<T>(curve, _t[i], _t[i + 2], _t[i + 1]);
    }
    _c[n] = _c[0];
  } else {
    for (int i = 0; i < n; i++) {
      _c[i] = new PSubCurve<T>(curve, _t[i], _t[i + 2], _t[i + 1]);
    }
  }
}

template <typename T>
void BlendingSpline<T>::eval(T t, int d, bool left) const {
  this->_p.setDim(d + 1);

  int i = _getIndex(t);

  auto b = B(t, i);

  auto c_m = c(i - 1, t);
  auto c_i = c(i, t);

  this->_p[0] = c_m + b[1] * (c_i - c_m);
}

template <typename T> int BlendingSpline<T>::_getIndex(T t) const {
  if (t == this->getEndP()) {
    return _t.size() - 3;
  } else {
    return std::distance(_t.begin(),
                         std::upper_bound(_t.begin(), _t.end(), t)) -
           1;
  }
}

template <typename T> T BlendingSpline<T>::w(T t, int i, int d) const {
  return (t - _t[i]) / (_t[i + d] - _t[i]);
}

template <typename T> Vector<T, 2> BlendingSpline<T>::B(T t, int i) const {
  T w_1 = _blend(w(t, i, 1));
  return Vector<T, 2>(1 - w_1, w_1);
}

template <typename T> T BlendingSpline<T>::_blend(T w) const {
  return w - T(1.0 / (M_2PI)) * sin(M_2PI * w);
}

template <typename T> Point<T, 3> BlendingSpline<T>::c(int i, T t) const {
  T w_t = w(t, i, 2);
  T c_t = _c[i]->getParEnd() * w_t + _c[i]->getParStart() * (1 - w_t);

  return _c[i]->evaluateParent(c_t, 0)[0];
}

template <typename T> void BlendingSpline<T>::showControlCurves() {

  for (PCurve<T, 3> *curve : _c) {
    curve->setCollapsed(true);
    curve->sample(20, 0);
    curve->toggleDefaultVisualizer();
    curve->setColor(GMcolor::aliceBlue());
    this->insert(curve);
  }
}

template <typename T>
void BlendingSpline<T>::localSimulate(double dt){
    double r = 0.000001 + (1.0 - (1.0/(48 * M_PI))*(27*cos(_angle) - 3*cos(3*_angle))) * dt;
    for (auto &curve: _c){
        curve->rotate(Angle(r), Vector<T,3>(0, 0, 1));
        curve->translate(Vector<T,3>(0, r, 0));
        curve->scale(Point<T,3>(1, 1, 1) * (1 + 0.00001 * sin(_angle)));
    }
    _angle += r;
    if (_angle >= M_2PI){_angle -= M_2PI;}
    this->resample();
    this->setEditDone();
}
} // namespace Custom
