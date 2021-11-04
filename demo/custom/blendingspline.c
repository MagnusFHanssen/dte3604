#include "blendingspline.h"

namespace Custom {

using namespace GMlib;

template <typename T>
BlendingSpline<T>::BlendingSpline(PCurve<T, 3> *copy, int n)
    : PCurve<T, 3>(20, 0, 0), _t_start{copy->getParStart()}, _t_end{copy->getParEnd()}, _closed{copy->isClosed()} {

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

  T frac = (end - start)/(n);

  for (int i = 1; i < n; i++) {
    _t.push_back(start + i * frac);
  }

  _t.push_back(end);
  _t.push_back(end);

  if (_closed) {
    _t[0] = start - frac;
    _t[n + 1] = end + frac;
  }
}

template <typename T>
void BlendingSpline<T>::_makeControlCurves(PCurve<T, 3> *curve) {
  const int n = _t.size() - 2;
  _c.resize(n);
  if (_closed) { // If the curve is closed, the ends are the same
    _c[0] = new PSubCurve<T>(curve, _t[n - 1], _t[2], _t[1]);
    for (int i = 1; i < n - 3; i++) {
      _c[i] = new PSubCurve<T>(curve, _t[i], _t[i + 2], _t[i + 1]);
    }
    _c[n - 3] = _c[0];
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

  this->_p[0] = b[0] * c(i - 1, t) + b[1] * c(i, t);
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

  return _c[i]->evaluateParent(c_t)[0];
}

template <typename T>
void BlendingSpline<T>::showControlCurves(){
    for (PCurve<T,3>* curve: _c){
        curve->sample(20, 0);
        curve->toggleDefaultVisualizer();
        curve->setColor(GMcolor::aliceBlue());
        this->insert(curve);
    }
}
} // namespace Custom
