#ifndef HOUGH_TRANSFORM_H
#define HOUGH_TRANSFORM_H

#include "utils.h"
#include <vector>

template <typename R_T, typename THETA_T>
struct HoughTransformer2d {
private:
  using uint32_t = uint;
  using traits = math_traits<R_T, THETA_T>;
public:
  HoughTransformer2d(R_T rStep, THETA_T thetaStep) : rStep(rStep), thetaStep(thetaStep) {}

  HoughSpace<R_T, THETA_T> transform(const std::vector<Point<R_T>> &points) const {
    auto sizeTheta = static_cast<size_t>(static_cast<THETA_T>(2) * traits::pi() / thetaStep) + 1;
    R_T maxR = 0;
    for (const auto &p : points) {
      R_T sqr = p.x * p.x + p.y * p.y;
      if (sqr == 0) continue;
      R_T arg = p.x / traits::sqrt(sqr);
      THETA_T t = (traits::pi() / static_cast<THETA_T>(2)) - traits::arcsin(arg);
      maxR = std::max(getR(p, t), maxR);
    }
    size_t sizeR = static_cast<size_t>(maxR / rStep) + 1;
    HoughSpace<R_T, THETA_T> space(rStep, thetaStep, sizeR, sizeTheta);
    for (const auto &p : points) {
      for (size_t i = 0; i != sizeTheta; ++i) {
        THETA_T t = i * thetaStep;
        R_T r = getR(p, t);
        space.update(r, i);
      }
    }
    return space;
  }

private:

  const R_T rStep;
  const THETA_T thetaStep;
};

#endif // HOUGH_TRANSFORM_H
