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

# compile each .cpp into an exe
FILE(GLOB SRCS *.cpp)
FOREACH(SFILE ${SRCS})
    GET_FILENAME_COMPONENT(PROG_NAME ${SFILE} NAME_WE)
    ADD_EXECUTABLE(${PROG_NAME} ${SFILE})
    TARGET_LINK_LIBRARIES(${PROG_NAME} ${GMSH_LIBRARIES})    
ENDFOREACH()
