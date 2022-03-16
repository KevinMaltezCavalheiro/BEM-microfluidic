PROJECT(MATH0471 CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

# ------------------------------------------------------------------------------
# Find libraries and setup compiler
# ------------------------------------------------------------------------------

# build type is "" by default in Linux
IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF()

# enable C++11
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

IF(APPLE)
    # on macOS, do not give priority to frameworks/apps
    SET(CMAKE_FIND_APPBUNDLE LAST)
    SET(CMAKE_FIND_FRAMEWORK LAST)
    SET(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architectures for Mac OS X" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp -I/opt/homebrew/include")
    list(APPEND EXTERNAL_LIBRARIES "-L/opt/homebrew/lib -lomp")
ENDIF()

# find gmsh-sdk
# gmsh.h
FIND_PATH(GMSH_INCLUDE_DIRS NAMES "gmsh.h")
MESSAGE(STATUS "GMSH_INCLUDE_DIRS=" ${GMSH_INCLUDE_DIRS})
if(NOT GMSH_INCLUDE_DIRS)
    MESSAGE(FATAL_ERROR "gmsh.h not found!")
ENDIF()
INCLUDE_DIRECTORIES(${GMSH_INCLUDE_DIRS})

# libgmsh.so
FIND_LIBRARY(GMSH_LIBRARIES gmsh)
MESSAGE(STATUS "GMSH_LIBRARIES=" ${GMSH_LIBRARIES})
IF(NOT GMSH_LIBRARIES)
    MESSAGE(FATAL_ERROR "gmsh library not found!")
ENDIF()

# find Eigen
find_path(EIGEN_INCLUDE_DIRS "Eigen/Dense" 
          PATHS "${PROJECT_SOURCE_DIR}/lib/eigen" "/usr/include/eigen3")
MESSAGE(STATUS "EIGEN_INCLUDE_DIRS=" ${EIGEN_INCLUDE_DIRS})
IF(NOT EIGEN_INCLUDE_DIRS)
    MESSAGE(FATAL_ERROR "Eigen include dir not found!")
ENDIF()
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})

# ------------------------------------------------------------------------------

FILE(GLOB SRCS *.cpp)
ADD_EXECUTABLE(solver ${SRCS})
TARGET_LINK_LIBRARIES(solver ${GMSH_LIBRARIES})
