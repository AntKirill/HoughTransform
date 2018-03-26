#ifndef HOUGH_TRANSFORM_H
#define HOUGH_TRANSFORM_H

#include "utils.h"
#include <vector>

template <typename R_T, typename THETA_T>
struct HoughSpace;

template <typename R_T, typename THETA_T>
struct HoughTransformer2d {
private:
  using uint32_t = uint;
  using traits = math_traits<R_T, THETA_T>;
public:
  HoughTransformer2d(R_T rStep, THETA_T thetaStep) : rStep(rStep), thetaStep(thetaStep) {}

  HoughSpace<R_T, THETA_T> transform(const std::vector<Point<R_T>> &points) const {
    auto sizeTheta = static_cast<size_t>(static_cast<THETA_T>(2) * traits::pi() / thetaStep) + 2;
    R_T maxR = 0;
    for (const auto &p : points) {
      R_T sqr = p.x * p.x + p.y * p.y;
      maxR = std::max(sqr, maxR);
    }
    maxR = traits::sqrt(maxR);
    size_t sizeR = static_cast<size_t>(maxR / rStep) + 10;
    HoughSpace<R_T, THETA_T> space(rStep, thetaStep, sizeR, sizeTheta);
    for (const auto &p : points) {
      for (size_t i = 0; i != sizeTheta; ++i) {
        R_T r = space.getR(p, i);
        space.update(r, i);
      }
    }
    return space;
  }

private:

  const R_T rStep;
  const THETA_T thetaStep;
};

template <typename R_T, typename THETA_T>
struct HoughSpace {
  const R_T rStep;
  const THETA_T thetaStep;

  using space_t = std::vector<std::vector<uint32_t>>;

  uint32_t get(R_T r, THETA_T theta) const {
    bool ok = true;
    Cell cell = getCell(r, theta, ok);
    if (!ok) return 0;
    return space[cell.rTimes][cell.thetaTimes];
  }
  space_t getSpace() const { return space; }

  std::vector<Line<R_T, THETA_T>> getLines(uint32_t amount) const {
    std::multiset<Node> lines;
    for (size_t rt = 0; rt < space.size(); ++rt) {
      for (size_t thetat = 0; thetat < space[0].size(); ++thetat) {
//        if ((rt == 0) && (thetat * thetaStep) >= math_traits<R_T, THETA_T>::pi()) continue;
        R_T liner = rHead[rt];
        THETA_T linetheta = thetaHead[thetat];
        if (space[rt][thetat] > 1) lines.emplace(space[rt][thetat], Line<R_T, THETA_T>(liner,
              linetheta));
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
    bool ok = true;
    Cell cellLine = getCell(line.r, line.theta, ok);
    assert(ok == true);
    size_t cell_theta = getCellComponent(line.theta, thetaStep, ok);
    R_T r = this->getR(p, cell_theta);
    Cell cellLine1 = getCell(r, line.theta, ok);
    if (!ok) return false;
    return cellLine == cellLine1;
  }

private:

  R_T getR(const Point<R_T> &p, size_t cell_theta) const {
    return ::getR(p, thetaHead[cell_theta]);
  }

  struct Cell {
    size_t rTimes, thetaTimes;
    Cell() : rTimes(0), thetaTimes(0) {}
    Cell(size_t i, size_t j) : rTimes(i), thetaTimes(j) {}
    Cell(Cell &&rhs) : rTimes(std::move(rhs.rTimes)), thetaTimes(std::move(rhs.thetaTimes)) {}
    Cell(const Cell &rhs) : rTimes(rhs.rTimes), thetaTimes(rhs.thetaTimes) {}
    bool operator==(const Cell &rhs) {
      return ((rTimes == rhs.rTimes) &&
              (thetaTimes == rhs.thetaTimes));
    }
  };

  template <typename T>
  size_t getCellComponent(T x, const T step, bool &ok) const {
    if (x < 0) {
      ok = false;
      return 0;
    }
    ok = true;
    auto cell_x = static_cast<size_t>(x / step);
    T newX = std::round(static_cast<T>(cell_x)) * step;
    T delta = x - newX;
    if ((0 < delta) && (delta < step)) return cell_x;
    if ((delta < 0) && (delta > -step)) return cell_x;
    if (delta > step) return cell_x + 1;
    if (delta < -step) return cell_x - 1;
    return cell_x;
  }

  Cell getCell(R_T r, THETA_T theta, bool &ok) const {
    size_t cell_r = getCellComponent(r, rStep, ok);
    if (!ok) return Cell();
    size_t cell_theta = getCellComponent(theta, thetaStep, ok);
    if (!ok) return Cell();
    return Cell(cell_r, cell_theta);
  }

  struct Node {
    uint32_t count;
    Line<R_T, THETA_T> line;
    bool operator<(const Node &rhs) const { return count < rhs.count; }
    Node(uint32_t count, Line<R_T, THETA_T> line) : count(count), line(line) {}
  };

  std::vector<std::vector<uint32_t>> space;
  std::vector<THETA_T> thetaHead;
  std::vector<R_T> rHead;

  HoughSpace(R_T rStep, THETA_T thetaStep, uint32_t rSize, uint32_t thetaSize) : rStep(rStep),
    thetaStep(thetaStep), space(rSize + 2, std::vector<uint32_t>(thetaSize + 2, 0)),
    thetaHead(thetaSize + 2), rHead(rSize + 2) {
    THETA_T thetaStep2 = thetaStep / static_cast<THETA_T>(2);
    R_T rStep2 = rStep / static_cast<R_T>(2);
    for (size_t i = 0; i < thetaSize; ++i) {
      thetaHead[i] = i * thetaStep + thetaStep2;
    }
    for (size_t i = 0; i < rSize; ++i) {
      rHead[i] = i * rStep + rStep2;
    }
  }

  friend HoughSpace<R_T, THETA_T> HoughTransformer2d<R_T, THETA_T>::transform(
    const std::vector<Point<R_T>> &) const;

  void update(R_T r, THETA_T theta) {
    bool ok = true;
    size_t cell_theta = getCellComponent(theta, thetaStep, ok);
    if (!ok) return;
    update(r, cell_theta);
  }

  void checkDist(size_t rt, size_t thetat) {
    if (rt >= space.size())
      std::cerr << rt << " >= " << space.size() << std::endl;
    if (thetat >= space[0].size())
      std::cerr << thetat << " >= " << space[0].size() << std::endl;
  }

  void update(size_t rt, size_t thetat) {
    #ifndef NDEBUG
    checkDist(rt, thetat);
    #endif
    ++space[rt][thetat];
  }

  void update(R_T r, size_t thetat) {
    bool ok = true;
    size_t cell_r = getCellComponent(r, rStep, ok);
    if (!ok) return;
    update(cell_r, thetat);
  }
};

#endif // HOUGH_TRANSFORM_H
