#include "csv.h"
#include <string>
#include <vector>

void get_parameter(std::string filename, std::vector<std::vector<double>>& data,
                   const int dim) {

  data.erase(data.begin());
  if (dim == 1) {
    io::CSVReader<2> in(filename);
    double v1, v2;
    while (in.read_row(v1, v2)) {
      data.push_back({v1, v2});
    }
  } else if (dim == 2) {
    io::CSVReader<3> in(filename);
    double v1, v2, v3;
    while (in.read_row(v1, v2, v3)) {
      data.push_back({v1, v2, v3});
    }
  } else if (dim == 3) {
    io::CSVReader<4> in(filename);
    double v1, v2, v3, v4;
    while (in.read_row(v1, v2, v3, v4)) {
      data.push_back({v1, v2, v3, v4});
    }
  }
}
