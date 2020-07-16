#include "boundaryvalues.h"
#include "clock.h"
#include "externalfunc.h"
#include "globalvariables.h"
#include "initialvalues.h"
#include "sourceterm.h"

#include "geothermal.h"

int main() {
  double seconds;
  get_parameter("inputfiles/parameters.csv", EquationData::g_perm_list, 1);

  if (EquationData::time_unit == 'd') {
    seconds = 86400;
  } else if (EquationData::time_unit == 'h') {
    seconds = 3600;
  }
  for (unsigned int i = 0; i < EquationData::g_time_sequence.size(); ++i) {
    EquationData::g_time_sequence[i] *= seconds;
  }

  try {
    using namespace dealii;
    CoupledTH<3> coupled_TH_solver(1);
    coupled_TH_solver.run();
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
