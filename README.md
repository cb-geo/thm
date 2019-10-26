# High-Performance Thermo-Hydro-Mechanical Code (CB-Geo thm)
> [CB-Geo Computational Geomechanics Research Group](https://www.cb-geo.com)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/thm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/thm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://thm.cb-geo.com/)
[![CircleCI](https://circleci.com/gh/cb-geo/thm.svg?style=svg)](https://circleci.com/gh/cb-geo/thm)
[![codecov](https://codecov.io/gh/cb-geo/thm/branch/develop/graph/badge.svg)](https://codecov.io/gh/cb-geo/thm)
[![](https://img.shields.io/github/issues-raw/cb-geo/thm.svg)](https://github.com/cb-geo/thm/issues)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/cb-geo/thm/projects/)

## Documentation

Please refer to [CB-Geo thm Documentation](https://cb-geo.github.io/thm-doc) for information on compiling, and running the code. The documentation also include the thm theory.

## Install dependencies

* Docker image for CB-Geo thm code [https://hub.docker.com/r/cbgeo/thm](https://hub.docker.com/r/cbgeo/thm)

* Instructions for running thm docker container: [https://github.com/cb-geo/docker-thm/blob/master/README.md](https://github.com/cb-geo/thm-container/blob/master/README.md).

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Boost](http://www.boost.org/)
* [Dealii](https://dealii.org)
* [Eigen](http://eigen.tuxfamily.org/)
* [Intel TBB](https://www.threadingbuildingblocks.org/)
* [HDF5](https://support.hdfgroup.org/HDF5/)

#### Optional

## Compile
> See https://thm-doc.cb-geo.com/ for more detailed instructions. 

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++ ..`.

1. Run `make clean && make -jN` (where N is the number of cores).
