#include "boundaryvalues.h"
#include "clock.h"
#include "externalfunc.h"
#include "globalvariables.h"
#include "initialvalues.h"
#include "sourceterm.h"

#include "geothermal.h"

int main(int argc, char** argv) {
  double seconds;
  // get_parameter("inputfiles/parameters.csv", EquationData::g_perm_list, 1);

  try {
    using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

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
