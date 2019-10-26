project(thm LANGUAGES CXX)

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

# CMake Modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})



# Include directories
include_directories(BEFORE
  ${thm_SOURCE_DIR}/include/
  ${thm_SOURCE_DIR}/external/
)

# thm executable
#SET(thm_src
#  ${thm_SOURCE_DIR}/src/thm.cc
#)

#add_library(lthm SHARED ${thm_src})

add_executable(thm ${thm_SOURCE_DIR}/src/main.cc)
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
