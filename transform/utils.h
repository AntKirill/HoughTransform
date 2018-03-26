#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <stdint.h>
#include "hough_transform.h"
#include <set>
#include <cmath>
#include <cassert>

template <typename R_T, typename THETA_T>
struct math_traits {
  static constexpr THETA_T pi() { return static_cast<THETA_T>(3.14159265359); }
  static R_T sin(const THETA_T &v) { return static_cast<R_T>(std::sin(v)); }
  static THETA_T arcsin(const R_T &v) { return static_cast<THETA_T>(std::asin(v)); }
  static R_T sqrt(const R_T &v) { return static_cast<R_T>(std::sqrt(v)); }
  static R_T cos(const THETA_T &v) { return static_cast<R_T>(std::cos(v)); }
};

template <typename T>
struct Point {
  const T x, y;
  Point(T x, T y) : x(x), y(y) {}
  Point(const Point<T> &rhs) : x(rhs.x), y(rhs.y) {}
};

template <typename R_T, typename THETA_T>
struct Line {
  const R_T r;
  const THETA_T theta;
  Line(R_T r, THETA_T theta) : r(r), theta(theta) {}
};

/*
 * Function to get distance between (0, 0) and line
 * @param p - any point on line
 * @param theta - angle between OX and normal to line from point (0, 0)
 * @return distance between point (0, 0) and line. If return value is negativ
 *  then such line does not exist
*/
template <typename R_T, typename THETA_T>
R_T getR(const Point<R_T> &p, const THETA_T &theta) {
  using traits = math_traits<R_T, THETA_T>;
  return p.x * traits::cos(theta) + p.y * traits::sin(theta);
}

#endif // UTILS_H




