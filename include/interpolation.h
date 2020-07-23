
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/tria.h>

#include <deal.II/base/function_lib.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <fstream>
#include <iostream>

using namespace dealii;

// Class Interpolation will create an instance of Interpolation used in class
// CoupedTH. It is constructed through tensor products from vectors. The number
// of vectors consisting of the instance will depend on the dim. Its dim will
// depend on the format of input data files. For example current input file is
// parameters_for_interpolation.txt and it provides the coordinates in three
// dimensions(x,y,z). So, the dim should be 3. In the example, the instance of
// class Interpolation is constructed through 3 vectors in x, y and z directions
// individually. Lengths of vectors in x,y and z directions is needed to
// construct the instance.
template <int dim>
class Interpolation {
 public:
  // Constructors for the case dim=1, dim2 and dim3
  Interpolation(int data_size_side1, int data_size_side2, int data_size_side3,
                std::string filename);
  Interpolation(int data_size_side1, int data_size_side2, std::string filename);
  Interpolation(int data_size_side1, std::string filename);
  // get interpolation value at specific points
  // dim=1
  double value(const double x);
  // dim=2
  double value(const double x, const double z);
  // dim=3
  double value(const double x, const double y, const double z);

 public:
  int dimension;
  const Functions::InterpolatedTensorProductGridData<dim> interpolation_data;
  int dimension_data;
  static std::array<std::vector<double>, 3> coordinate_values(
      int data_size_side1, int data_size_side2, int data_size_side3,
      std::string filename);
  static std::array<std::vector<double>, 2> coordinate_values(
      int data_size_side1, int data_size_side2, std::string filename);
  static std::array<std::vector<double>, 1> coordinate_values(
      int data_size_side1, std::string filename);
  static std::vector<double> get_data(int data_size_side1, int data_size_side2,
                                      int data_size_side3,
                                      std::string filename);
  static std::vector<double> get_data(int data_size_side1, int data_size_side2,
                                      std::string filename);
  static std::vector<double> get_data(int data_size_side1,
                                      std::string filename);
  static Table<3, double> table_get_data(int data_size_side1,
                                         int data_size_side2,
                                         int data_size_side3,
                                         std::string filename);
  static Table<2, double> table_get_data(int data_size_side1,
                                         int data_size_side2,
                                         std::string filename);
  static Table<1, double> table_get_data(int data_size_side1,
                                         std::string filename);

 public:
  int size_side1;
  int size_side2;
  int size_side3;
  std::string input_data_filename;
};
template <int dim>
Interpolation<dim>::Interpolation(int data_size_side1, int data_size_side2,
                                  int data_size_side3, std::string filename)
    : interpolation_data(coordinate_values(data_size_side1, data_size_side2,
                                           data_size_side3, filename),
                         table_get_data(data_size_side1, data_size_side2,
                                        data_size_side3, filename)),
      size_side1(data_size_side1),
      size_side2(data_size_side2),
      size_side3(data_size_side3),
      input_data_filename(filename) {}

template <int dim>
Interpolation<dim>::Interpolation(int data_size_side1, int data_size_side2,
                                  std::string filename)
    : interpolation_data(
          coordinate_values(data_size_side1, data_size_side2, filename),
          table_get_data(data_size_side1, data_size_side2, filename)),
      size_side1(data_size_side1),
      size_side2(data_size_side2),
      size_side3(0),
      input_data_filename(filename) {}

template <int dim>
Interpolation<dim>::Interpolation(int data_size_side1, std::string filename)
    : interpolation_data(coordinate_values(data_size_side1, filename),
                         table_get_data(data_size_side1, filename)),
      size_side1(data_size_side1),
      size_side2(0),
      size_side3(0),
      input_data_filename(filename) {}

template <int dim>
std::array<std::vector<double>, 3> Interpolation<dim>::coordinate_values(
    int data_size_side1, int data_size_side2, int data_size_side3,
    std::string filename) {
  std::array<std::vector<double>, dim> coordinate_valuesoutput;
  std::vector<double> data_x;
  std::vector<double> data_y;
  std::vector<double> data_z;
  // create a stream where we read from data in the folder inputfiles
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  for (int line = 1;
       line <= data_size_side1 * data_size_side2 * data_size_side3; ++line) {
    double x, y, z, permeability;
    in >> x >> y >> z >> permeability;
    if (line <= data_size_side3) {
      data_z.push_back(z);  // push coordinate
    }
    if (line % data_size_side3 == 1 &&
        line <= data_size_side2 * data_size_side3) {
      data_y.push_back(y);  // push coordinate
    }
    if ((line % (data_size_side2 * data_size_side3)) == 1) {
      data_x.push_back(x);  // push coordinate
    }
  }
  coordinate_valuesoutput[0] = data_x;
  coordinate_valuesoutput[1] = data_y;
  coordinate_valuesoutput[2] = data_z;
  return coordinate_valuesoutput;
}

template <int dim>
std::array<std::vector<double>, 2> Interpolation<dim>::coordinate_values(
    int data_size_side1, int data_size_side2, std::string filename) {
  std::array<std::vector<double>, dim> coordinate_valuesoutput;
  std::vector<double> data_x;
  std::vector<double> data_z;
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  for (int line = 1; line <= data_size_side1 * data_size_side2; ++line) {
    double x, z, permeability;
    in >> x >> z >> permeability;
    if (line <= data_size_side2) {
      data_z.push_back(z);  // push coordinate
    }
    if (line % data_size_side1 == 1) {
      data_x.push_back(x);  // push coordinate
    }
  }
  coordinate_valuesoutput[0] = data_x;
  coordinate_valuesoutput[1] = data_z;
  return coordinate_valuesoutput;
}

template <int dim>
std::array<std::vector<double>, 1> Interpolation<dim>::coordinate_values(
    int data_size_side1, std::string filename) {
  std::array<std::vector<double>, dim> coordinate_valuesoutput;
  std::vector<double> data_z;
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  for (int line = 1; line <= data_size_side1; ++line) {
    double x, z, permeability;
    in >> z >> permeability;
    if (line <= data_size_side1) {
      data_z.push_back(z); // push coordinate
    }
  }
  coordinate_valuesoutput[0] = data_z;
  return coordinate_valuesoutput;
}

template <int dim>
std::vector<double> Interpolation<dim>::get_data(int data_size_side1,
                                                 int data_size_side2,
                                                 int data_size_side3,
                                                 std::string filename) {
  std::vector<double> data;
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  if (dim == 1) {
    for (int line = 0; line < data_size_side1; ++line) {
      double z, permeability;
      in >> z >> permeability;
      data.push_back(permeability); // push parameter
    }
    return data;
  } else if (dim == 2) {
    for (int line = 0; line < data_size_side1 * data_size_side2; ++line) {
      double x, z, permeability;
      in >> x >> z >> permeability;
      data.push_back(permeability);  // push parameter
    }
    return data;
  } else if (dim == 3) {
    for (int line = 0; line < data_size_side1 * data_size_side2; ++line) {
      double x, y, z, permeability;
      in >> x >> y >> z >> permeability;
      data.push_back(permeability);  // push parameter
    }
    return data;
  }
}


// get value of parameter
template <int dim>
std::vector<double> Interpolation<dim>::get_data(int data_size_side1,
                                                 int data_size_side2,
                                                 std::string filename) {
  std::vector<double> data;
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  for (int line = 0; line < data_size_side1 * data_size_side2; ++line) {
    double x, z, permeability;
    in >> x >> z >> permeability;
    data.push_back(permeability);  // push parameters
  }
  return data;
}

template <int dim>
std::vector<double> Interpolation<dim>::get_data(int data_size_side1,
                                                 std::string filename) {
  std::vector<double> data;
  boost::iostreams::filtering_istream in;
  in.push(boost::iostreams::file_source(filename));
  for (int line = 0; line < data_size_side1; ++line) {
    double x, z, permeability;
    in >> z >> permeability;
    data.push_back(permeability); // push parameters
  }
  return data;
}

template <int dim>
double Interpolation<dim>::value(const double x, const double y,
                                 const double z) {
  return interpolation_data.value(Point<3>(x, y, z));
}

template <int dim>
double Interpolation<dim>::value(const double x, const double z) {
  return interpolation_data.value(Point<2>(x, z));
}

template <int dim>
double Interpolation<dim>::value(const double z) {
  return interpolation_data.value(Point<1>(z));
}

template <int dim>
Table<3, double> Interpolation<dim>::table_get_data(int data_size_side1,
                                                    int data_size_side2,
                                                    int data_size_side3,
                                                    std::string filename) {
  return Table<3, double>(
      data_size_side1, data_size_side2, data_size_side3,
      get_data(data_size_side1, data_size_side2, data_size_side3, filename)
          .begin());
}

template <int dim>
Table<2, double> Interpolation<dim>::table_get_data(int data_size_side1,
                                                    int data_size_side2,
                                                    std::string filename) {
  return Table<2, double>(
      data_size_side1, data_size_side2,
      get_data(data_size_side1, data_size_side2, filename).begin());
}

template <int dim>
Table<1, double> Interpolation<dim>::table_get_data(int data_size_side1,
                                                    std::string filename) {
  return Table<1, double>(data_size_side1,
                          get_data(data_size_side1, filename).begin());
}
