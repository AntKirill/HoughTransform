#include <gtest.h>
#include "hough_transform.h"
#include "utils.h"

namespace {
  template <typename T1, typename T2>
  inline void checkEachLine(const std::vector<Line<T1, T2>> &lines, const std::vector<Point<T1>> &points,
                     const HoughSpace<T1, T2> &hs) {
    for (const auto &line : lines) {
      size_t cnt = 0;
      for (const auto &p : points) {
        if (hs.isOnLine(line, p)) ++cnt;
      }
      EXPECT_EQ(cnt, hs.get(line.r, line.theta));
    }
  }
}

TEST(small, test1) {
  std::vector<Point<double>> points = {{0, 0}, {1, 1}};
  HoughTransformer2d<double, double> transformer(0.1, 0.1);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
}

TEST(small, test2) {
  std::vector<Point<double>> points = {{0, 0}, {-1, -1}, {2, 2}};
  HoughTransformer2d<double, double> transformer(0.1, 0.001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_LE(2, hs.get(lines[0].r, lines[0].theta));
}

TEST(small, test3) {
  std::vector<Point<double>> points = {{0, 0.6666}, {-1, 0}, {2, 2}};
  HoughTransformer2d<double, double> transformer(0.1, 0.0001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_LE(2, hs.get(lines[0].r, lines[0].theta));
}


