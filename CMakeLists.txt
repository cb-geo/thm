#project(thm LANGUAGES CXX)
# Set the name of the project and target:
SET(TARGET "thm")

# Require C++14-compliant compiler; only available for CMake v3.10 and up
set(CMAKE_CXX_STANDARD 14)

cmake_minimum_required(VERSION 3.10)

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)

# General compile settings
IF (NOT CMAKE_BUILD_TYPE)
 SET(CMAKE_BUILD_TYPE "Debug")
 #SET(CMAKE_BUILD_TYPE "Release")
ENDIF (NOT CMAKE_BUILD_TYPE)

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# CMake seems to have no way to enable/disable testing per subproject,
# so we provide an option similar to BUILD_TESTING, but just for THM.
option(THM_BUILD_TESTING "enable testing for thm" ON)

# Add -O0 to remove optimizations when using gcc
# IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
  SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
# ENDIF(CMAKE_COMPILER_IS_GNUCC)

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Dealii
FIND_PACKAGE(deal.II 8.5.0 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
if(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

if (${deal.II_FOUND})
  include_directories(${DEAL_II_DIR}/include/)
endif()

#thm executable
SET(thm_src
  ./src/main.cc
  ./src/get_parameter.cc
  ./src/interpolation1d.cc
  ./src/linspace.cc
  ./src/inOrNot.cc
)

# aux_source_directory(./src/ thm_src)

# Include directories
include_directories(
  ./include/
  ./external/
  
)

SET(TARGET_SRC ${thm_src})

IF(NOT DEAL_II_WITH_PETSC OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
  MESSAGE(FATAL_ERROR "
Error! The deal.II library found at ${DEAL_II_PATH} was not configured with
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
One or all of these are OFF in your installation but are required for this tutorial step."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()


#add_library(lthm SHARED ${thm_src})

# add_executable(geothm ${thm_SOURCE_DIR}/src/main.cc)
#target_link_libraries(thm lthm)

# Unit test
#if(THM_BUILD_TESTING)
#  SET(test_src
#    ${thm_SOURCE_DIR}/tests/test_main.cc
#  )   
#  add_executable(thmtest ${test_src})
#  target_link_libraries(thmtest lthm)
#  add_test(NAME thmtest COMMAND $<TARGET_FILE:thmtest>)
#  enable_testing()
#endif()

# Coverage
# find_package(codecov)
# if(ENABLE_COVERAGE)
#  add_executable(thmtest_coverage ${thm_src} ${test_src})
#  add_coverage(thmtest_coverage)
#endif()
