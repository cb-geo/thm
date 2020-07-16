#include "csv.h"
#include <string>
#include <vector>

void get_parameter(std::string filename, std::vector<std::vector<double>>& data,
                   const int dim) {

  data.erase(data.begin());
  if (dim == 1) {
    io::CSVReader<2> in(filename);
    double variable1, variable2;
    while (in.read_row(variable1, variable2)) {
      data.push_back({variable1, variable2});
    }
  } else if (dim == 2) {
    io::CSVReader<3> in(filename);
    double variable1, variable2, variable3;
    while (in.read_row(variable1, variable2, variable3)) {
      data.push_back({variable1, variable2, variable3});
    }
  } else if (dim == 3) {
    io::CSVReader<4> in(filename);
    double variable1, variable2, variable3, variable4;
    while (in.read_row(variable1, variable2, variable3, variable4)) {
      data.push_back({variable1, variable2, variable3, variable4});
    }
  }
}
