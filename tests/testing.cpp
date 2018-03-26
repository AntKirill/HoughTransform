#include <gtest.h>
#include "hough_transform.h"
#include "utils.h"
#include <random>

namespace {
  template <typename T1, typename T2>
  inline void checkEachLine(const std::vector<Line<T1, T2>> &lines,
                            const std::vector<Point<T1>> &points,
                            const HoughSpace<T1, T2> &hs) {
    for (const auto &line : lines) {
      size_t cnt = 0;
      for (const auto &p : points) {
        if (hs.isOnLine(line, p)) ++cnt;
      }
      EXPECT_EQ(cnt, hs.get(line.r, line.theta));
    }
  }

  template<typename T>
  std::vector<Point<T>> generatePoints(int begAmount, int endAmount, const Point<T> &leftDown,
  const Point<T> &rightUp) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(begAmount, endAmount);
    std::uniform_real_distribution<T> disPointX(leftDown.x, rightUp.x);
    std::uniform_real_distribution<T> disPointY(leftDown.y, rightUp.y);
    int n = dis(gen);
    std::vector<Point<float>> points;
    for (int i = 0; i < n; ++i) {
      points.push_back({disPointX(gen), disPointY(gen)});
    }
    return points;
  }

}

TEST(twoPoints, test1) {
  std::vector<Point<double>> points = {{0, 0}, {1, 1}};
  HoughTransformer2d<double, double> transformer(0.1, 0.1);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
}

TEST(twoPoints, test2) {
  std::vector<Point<double>> points = {{0, 0}, {1, 1}};
  HoughTransformer2d<double, float> transformer(0.1, 0.1);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
}

TEST(threePoints, test1) {
  std::vector<Point<double>> points = {{0, 0}, {-1, -1}, {2, 2}};
  HoughTransformer2d<double, double> transformer(0.1, 0.001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_LE(2, hs.get(lines[0].r, lines[0].theta));
}

TEST(threePoints, test2) {
  std::vector<Point<double>> points = {{0, 0.66}, {-1, 0}, {2, 2}};
  HoughTransformer2d<double, double> transformer(0.01, 0.0001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_LE(2, hs.get(lines[0].r, lines[0].theta));
}

TEST(threePoints, test3) {
  std::vector<Point<double>> points = {{0, 0.6666666666}, {-1, 0}, {2, 2}};
  HoughTransformer2d<double, double> transformer(0.001, 0.0001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_EQ(3, hs.get(lines[0].r, lines[0].theta));
}

TEST(threePoints, test4) {
  std::vector<Point<float>> points = {{0, 0.6666666}, {-1, 0}, {2, 2}};
  HoughTransformer2d<float, long double> transformer(0.001, 0.0001);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000);
  checkEachLine(lines, points, hs);
  EXPECT_EQ(3, hs.get(lines[0].r, lines[0].theta));
}

TEST(manyPoints, rarePoints) {
  auto points = generatePoints<float>(10, 20, Point<float>(-100, -100), Point<float>(100, 100));
  HoughTransformer2d<float, float> transformer(0.01, 0.01);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000000000);
  checkEachLine(lines, points, hs);
}

TEST(manyPoints, frequentPoints) {
  auto points = generatePoints<float>(100, 110, Point<float>(-1, -1), Point<float>(1, 1));
  HoughTransformer2d<float, float> transformer(0.01, 0.01);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000000000);
  checkEachLine(lines, points, hs);
}

TEST(manyPoints, randomPointsDestribution) {
  int N = 10;
  while (N--) {
    auto points = generatePoints<float>(0, 1000, Point<float>(-10, -10), Point<float>(10, 10));
    HoughTransformer2d<float, float> transformer(0.01, 0.01);
    auto hs = transformer.transform(points);
    auto lines = hs.getLines(1000000000);
    checkEachLine(lines, points, hs);
  }
}

TEST(Points10000, randomPointsDestribution) {
  auto points = generatePoints<float>(10000, 10000,
                                      Point<float>(-1000, -1000), Point<float>(1000, 1000));
  HoughTransformer2d<float, float> transformer(0.01, 0.01);
  auto hs = transformer.transform(points);
  auto lines = hs.getLines(1000000000);
  checkEachLine(lines, points, hs);
}
