#include "geothermal.h"
int main() {
  try {
    using namespace dealii;
    using namespace Geothermal;
    HeatEquation<3> heat_equation_solver;
    heat_equation_solver.run();
  } catch (std::exception& exc) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  } catch (...) {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  return 0;

  // try
  //   {
  //     grid_input ();
  //   }
  // catch (std::exception &exc)
  //   {
  //     std::cerr << std::endl << std::endl
  //               << "----------------------------------------------------"
  //               << std::endl;
  //     std::cerr << "Exception on processing: " << std::endl
  //               << exc.what() << std::endl
  //               << "Aborting!" << std::endl
  //               << "----------------------------------------------------"
  //               << std::endl;

  //     return 1;
  //   }
  // catch (...)
  //   {
  //     std::cerr << std::endl << std::endl
  //               << "----------------------------------------------------"
  //               << std::endl;
  //     std::cerr << "Unknown exception!" << std::endl
  //               << "Aborting!" << std::endl
  //               << "----------------------------------------------------"
  //               << std::endl;
  //     return 1;
  //   }
}