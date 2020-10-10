#include "boundaryvalues.h"
#include "clock.h"
#include "externalfunc.h"
#include "globalvariables.h"
#include "initialvalues.h"
#include "sourceterm.h"
#include "geothermal.h"
// edit
//#include "interpolation.h"

int main(int argc, char** argv) {
  using namespace dealii;

  if (EquationData::time_unit == 'd') {
    EquationData::g_seconds = 86400;
  } else if (EquationData::time_unit == 'h') {
    EquationData::g_seconds = 3600;
  } else if (EquationData::time_unit == 's') {
    EquationData::g_seconds = 1;
  }

  for (unsigned int i = 0; i < EquationData::g_time_sequence.size(); ++i) {
    EquationData::g_time_sequence[i] *= EquationData::g_seconds;
  }
  // edit
  // std::cout << 1<< std::endl;
  // std::cout << data_interpolation.value(0.0, 0.0, 0.0) << std::endl;

  try {
    // edit
    // using namespace dealii;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

    // edit
    CoupledTH<3> coupled_TH_solver(1);
    coupled_TH_solver.run();
    // std::string filename2="inputfiles/parameters_for_interpolation.txt";
    // Interpolation<3> output_data2=Interpolation<3>(3,3,3,filename2);
    // std::cout<<output_data2.value(-1,-1,-1)<<std::endl;
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
