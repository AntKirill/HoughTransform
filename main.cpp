#include <iostream>
#include <fstream>
#include <hough_transform.h>
#include <utils.h>

using namespace std;

void printUsage(const string &str = "") {
  if (str != "") cout << str << endl;
  cout << "Usage:\nFirst argument - file name.";
}

/*
 * Gets amount one argument - file name.
 * File contains n - amount of points on first line, and next amount lines contains points,
 * and amount of lines on the last line.
 * Produce output on stdout with lines.
*/
int main(int argc, char *argv[]) {
  if (argc != 2) {
    printUsage();
    return 0;
  }
  ifstream fin;
  fin.open(argv[1]);
  if (!fin.good()) {
    string str = "file " + std::string(argv[1]) + " not found";
    printUsage(str);
    return 0;
  }
  int n;
  fin >> n;
  vector<Point<double>> points;
  for (int i = 0; i < n; ++i) {
    double x, y;
    fin >> x >> y;
    points.emplace_back(Point<double> {x, y});
  }

  HoughTransformer2d<double, double> transformer(0.1, 0.01);
  auto hs = transformer.transform(points);
  size_t m = 0;
  fin >> m;
  auto lines = hs.getLines(m);
  auto space = hs.getSpace();

//  cout << "HoughSpace: " << endl;
//  for (auto vi : space) {
//    for (auto i : vi) cout << i << " ";
//    cout << endl;
//  }
//  cout << endl;
  for (const auto &line : lines) {
    cout << line.r << " " << line.theta << endl;
    for (const auto &p : points)
      cout << hs.isOnLine(line, p) << endl;
  }
  cout << endl;

  cout << hs.get(lines[0].r, lines[0].theta) << endl;

  fin.close();
  return 0;
}
