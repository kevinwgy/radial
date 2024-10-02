cmake_minimum_required(VERSION 3.10)

# set the project name
project(radial VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

# compiler flags (turn on all the warnings)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#--------------------------------------------------------
# find packages
find_package(Eigen3 3.3 REQUIRED)
#--------------------------------------------------------

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
include_directories(${CMAKE_SOURCE_DIR}/MathTools)

# add the executable
add_executable(radial 
Main.cpp
MathTools/rbf_interp.cpp)
