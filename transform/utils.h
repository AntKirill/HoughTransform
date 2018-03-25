#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <stdint.h>
#include "hough_transform.h"
#include <set>
#include <cmath>

template <typename R_T, typename THETA_T>
struct HoughTransformer2d;

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

template <typename T>
struct myoptional {
  void set(const T &value) { this->value = value; }
  bool hasValue() const { return _hasValue; }
  T getValue() const {return value;}
  myoptional() : _hasValue(false) {}
  myoptional(T value) : _hasValue(true), value(value) {}
  myoptional(const myoptional<T> &rhs) { *this = rhs;}
  myoptional &operator=(const myoptional<T> &rhs) {
    _hasValue = rhs._hasValue;
    value = rhs.value;
    return *this;
  }
  myoptional &operator=(myoptional<T> &&rhs) {
    _hasValue = rhs._hasValue;
    value = std::move(rhs.value);
    return *this;
  }
private :
  bool _hasValue;
  T value;
};

template <typename R_T, typename THETA_T>
struct HoughSpace {
  const R_T rStep;
  const THETA_T thetaStep;

  using space_t = std::vector<std::vector<uint32_t>>;

  uint32_t get(R_T r, THETA_T theta) const {
    auto rt = static_cast<size_t>(r / rStep);
    auto thetat = static_cast<size_t>(theta / thetaStep);
    return space[rt][thetat];
  }
  space_t getSpace() const { return space; }

  std::vector<Line<R_T, THETA_T>> getLines(uint32_t amount) const {
    std::multiset<Node> lines;
    for (size_t rt = 0; rt < space.size(); ++rt) {
      for (size_t thetat = 0; thetat < space[0].size(); ++thetat) {
        //         if ((rt == 0) && (thetat * thetaStep) >= math_traits<R_T, THETA_T>::pi()) continue;
        if (space[rt][thetat] > 1) lines.emplace(space[rt][thetat], Line<R_T, THETA_T>(rt * rStep,
              thetat * thetaStep));
      }
    }
    std::vector<Line<R_T, THETA_T>> amountLines;
    uint32_t cnt = 0;
    for (auto it = lines.rbegin(); (it != lines.rend()) && (cnt != amount); ++it, ++cnt) {
      amountLines.emplace_back(std::move(it->line));
    }
    return amountLines;
  }

  bool isOnLine(const Line<R_T, THETA_T> &line, const Point<R_T> &p) const {
    Cell cellLine = getCell(line);
    R_T r = getR(p, line.theta);
    Cell cellLine1 = getCell(r, line.theta);
    return cellLine == cellLine1;
  }

private:

  struct Cell {
    myoptional<size_t> rTimes, thetaTimes;
    Cell() = default;
    Cell(int32_t i, int32_t j) {
      if (i >= 0) rTimes = myoptional<size_t>(static_cast<size_t>(i));
      if (j >= 0) thetaTimes = myoptional<size_t>(static_cast<size_t>(j));
    }
    Cell(Cell &&rhs) : rTimes(std::move(rhs.rTimes)), thetaTimes(std::move(rhs.thetaTimes)) {}
    Cell(const Cell &rhs) : rTimes(rhs.rTimes), thetaTimes(rhs.thetaTimes) {}
    bool operator==(const Cell &rhs) {
      return ((hasValues() && rhs.hasValues()) &&
              (rTimes.getValue() == rhs.rTimes.getValue()) &&
              (thetaTimes.getValue() == rhs.thetaTimes.getValue()));
    }
    bool hasValues() const { return rTimes.hasValue() && thetaTimes.hasValue(); }
  };

  Cell getCell(R_T r, THETA_T theta) const {
    if (r < 0) return Cell();
    auto rt = static_cast<size_t>(r / rStep);
    auto thetat = static_cast<size_t>(theta / thetaStep);
    return Cell(rt, thetat);
  }

  Cell getCell(R_T r, size_t theta) const {
    if (r < 0) return Cell();
    auto rt = static_cast<size_t>(r / rStep);
    return Cell(rt, theta);
  }

  Cell getCell(const Line<R_T, THETA_T> &line) const {
    return getCell(line.r, line.theta);
  }

  struct Node {
    uint32_t count;
    Line<R_T, THETA_T> line;
    bool operator<(const Node &rhs) const { return count < rhs.count; }
    Node(uint32_t count, Line<R_T, THETA_T> line) : count(count), line(line) {}
  };

  std::vector<std::vector<uint32_t>> space;

  HoughSpace(R_T rStep, THETA_T thetaStep, const space_t &space) : rStep(rStep), thetaStep(thetaStep),
    space(space) {}
  HoughSpace(R_T rStep, THETA_T thetaStep, uint32_t rSize, uint32_t thetaSize) : rStep(rStep),
    thetaStep(thetaStep), space(rSize, std::vector<uint32_t>(thetaSize, 0)) {}

  friend HoughSpace<R_T, THETA_T> HoughTransformer2d<R_T, THETA_T>::transform(
    const std::vector<Point<R_T>> &) const;

  void update(const R_T &r, const THETA_T &theta) {
    Cell cell = getCell(r, theta);
    if (!cell.hasValues()) return;
    update(cell.rTimes.getValue(), cell.thetaTimes.getValue());
  }

  void update(size_t rt, size_t thetat) {
    ++space[rt][thetat];
  }

  void update(const R_T &r, size_t thetat) {
   Cell cell = getCell(r, thetat);
   if (!cell.hasValues()) return;
   update(cell.rTimes.getValue(), cell.thetaTimes.getValue());
  }

};

#endif // UTILS_H




