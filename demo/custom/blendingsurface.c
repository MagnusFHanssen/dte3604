#include "blendingsurface.h"

namespace Custom {

using namespace GMlib;

template <typename T>
inline BlendingSurface<T>::BlendingSurface(PSurf<T, 3> *copy, int u, int v)
    : PSurf<T, 3>(), _startU{copy->getParStartU()}, _endU{copy->getParEndU()},
      _startV{copy->getParStartV()}, _endV{copy->getParEndV()},
      _closedU{copy->isClosedU()}, _closedV{copy->isClosedV()}, _surface{copy} {
  _fillKnotVector(_tU, _startU, _endU, u, _closedU);
  _fillKnotVector(_tV, _startV, _endV, v, _closedV);

  _makeLocalSurfaces(copy);
}

template <typename T> BlendingSurface<T>::~BlendingSurface() {}

// Getters
template <typename T> T BlendingSurface<T>::getStartPU() const {
  return _startU;
}

template <typename T> T BlendingSurface<T>::getEndPU() const { return _endU; }

template <typename T> T BlendingSurface<T>::getStartPV() const {
  return _startV;
}

template <typename T> T BlendingSurface<T>::getEndPV() const { return _endV; }

template <typename T> bool BlendingSurface<T>::isClosedU() const {
  return _closedU;
}

template <typename T> bool BlendingSurface<T>::isClosedV() const {
  return _closedV;
}

template <typename T>
void BlendingSurface<T>::_fillKnotVector(std::vector<T> &t, const T startP,
                                         const T endP, int steps, bool closed) {
  t = {startP, startP};

  int n = steps - (closed ? 0 : 1);

  T frac = (endP - startP) / n;

  for (int i = 1; i < n; i++) {
    t.push_back(startP + i * frac);
  }

  t.push_back(endP);
  t.push_back(endP);

  if (closed) {
    t[0] = startP - frac;
    t[n + 2] = endP + frac;
  }
}

template <typename T>
void BlendingSurface<T>::_makeLocalSurfaces(PSurf<T, 3> *surface) {
  _s.resize(_tU.size() - _closedU);

  for (int i = 0; i < _tU.size() - 2 - _closedU; i++) {
    _s[i].resize(_tV.size() - _closedV);
    for (int j = 0; j < _tV.size() - 2 - _closedV; j++) {
      _s[i][j] = new SubPatch<T>(surface, _tU[i + 1], _tV[i + 1]);
    }
    if (_closedV) {
      _s[i][_tV.size() - 3] = _s[i][0];
    }
  }
  if (_closedU) {
    _s[_tU.size() - 3] = _s[0];
  }
}

template <typename T> T BlendingSurface<T>::_getIndex(T t, bool is_u) const {
    std::vector<T> _t{(is_u ? _tU : _tV)};
    if (t == (is_u ? this->getEndPU() : this->getEndPV())) {
    return _t.size() - 3;
  } else {
    return std::distance(_t.begin(),
                         std::upper_bound(_t.begin(), _t.end(), t)) -
           1;
  }
}

template <typename T>
T BlendingSurface<T>::w(std::vector<T> const& knots, T const t, int i, int d) const {
  return (t - knots[i]) / (knots[i + d] - knots[i]);
}

template <typename T> Vector<T, 2> BlendingSurface<T>::_blend(T w) const {
    return {T(3) * w * w - T(2) * w * w * w, T(6) * w * (T(1) - w)};
}

template <typename T>
void BlendingSurface<T>::eval(T u, T v, int d1, int d2, bool lu,
                              bool lv) const {
    this->_p.setDim(d1 + 1, d2 + 1);
    int i_u = _getIndex(u, true);
    int i_v = _getIndex(v, false);

    HqMatrix<T,3> A_[2][2] = {
        {_s[i_u - 1][i_v - 1]->getMatrix(), _s[i_u - 1][i_v]->getMatrix()},
        {_s[i_u][i_v - 1]->getMatrix(), _s[i_u][i_v]->getMatrix()}
    };

    HqMatrix<T,3> A_arr[3] = {
        A_[1][0] - A_[0][0],
        A_[0][1] - A_[0][0],
        A_[1][1] + A_[0][0] - A_[0][1] - A_[1][0]
    };

    Vector<T,2> b_u = _blend(w(_tU, u, i_u, 1));
    Vector<T,2> b_v = _blend(w(_tV, v, i_v, 1));

    HqMatrix<T,3> A = A_[0][0] + A_arr[0] * b_u[0] + A_arr[1] * b_u[1] + A_arr[2] * (b_u[0] * b_v[0]);
    HqMatrix<T,3> A_u = A_arr[0] * b_u[1] + A_arr[2] * (b_u[1] * b_v[0]);
    HqMatrix<T,3> A_v = A_arr[1] * b_v[1] + A_arr[2] * (b_u[0] * b_v[1]);

    DMatrix<Vector<T,3>> S = _surface->evaluate(u, v, d1, d2);


    Point<T,3> p = S[0][0];
    Vector<T,3> S_u = S[1][0];
    Vector<T,3> S_v = S[0][1];

    this->_p[0][0] = A * p;
    this->_p[1][0] = A_u * p + A * S_u;
    this->_p[0][1] = A_v * p + A * S_v;

}

template <typename T>
void BlendingSurface<T>::showLocalSurfaces(){
    for (std::vector<SubPatch<T>*> vec : _s){
        for (SubPatch<T>* surf : vec){
            surf->setCollapsed(true);
            this->insert(surf);
        }
    }
}

} // namespace Custom
