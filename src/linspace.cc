#include <vector>

std::vector<double> linspace(double a, double b, int N) {
  double h = (b - a) / static_cast<double>(N - 1);
  std::vector<double> xs(N);
  typename std::vector<double>::iterator x;
  double val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = val;
  return xs;
}