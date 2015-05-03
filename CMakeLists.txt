# CMakeLists.txt

cmake_minimum_required(VERSION 2.6)
project(smartplus)

# The version number
set(VERSION_MAJOR 0)
set(VERSION_MINOR 9)
set(VERSION_PATCH 0)

if(VERSION_PATCH MATCHES "0")
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}")
else()
	set(VERSION_NUMBER "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
endif()

message(STATUS "Configuring ${PROJECT_NAME} ${VERSION_NUMBER}")

# We default to debugging mode for developers.
option(DEBUG "Compile with debugging information" OFF)

# configure a header file to pass some CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/configuration.hpp.in"
  "${PROJECT_BINARY_DIR}/configuration.hpp"
)

#Set executable files and library files
set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib/${CMAKE_BUILD_TYPE})

#Inclusion of Armadillo Boost
include(FindArmadillo)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# OpenMP
include(FindOpenMP)
#find_package(OpenMP)

# Build type
if(NOT CMAKE_BUILD_TYPE)  # Debug by default
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel"
        FORCE)
endif()

message(STATUS "Build type ${CMAKE_BUILD_TYPE}")

# Set the CFLAGS and CXXFLAGS depending on the options the user specified.
# Only GCC-like compilers support -Wextra, and other compilers give tons of
# output for -Wall, so only -Wall and -Wextra on GCC.
if(CMAKE_COMPILER_IS_GNUCC OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
endif(CMAKE_COMPILER_IS_GNUCC OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

# Debugging CFLAGS.  Turn optimizations off; turn debugging symbols on.
if(DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O0")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3")
endif(DEBUG)

if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


#Inclusion of public headers
include_directories(include)

# add the binary tree to the search path for include files
include_directories("${PROJECT_BINARY_DIR}")

#Command file to get all the files in the src/ and include/
file(GLOB_RECURSE source_files src/* include/*)
list(REMOVE_ITEM source_files "Control/solver.cpp")

#Add the files to the lib
add_library(smartplus SHARED ${source_files})
#link against armadillo
target_link_libraries(smartplus ${ARMADILLO_LIBRARIES})

#Add the solver executable
add_executable(solver src/Control/solver.cpp ${source_files})

#Link the solver executable with armadillo
target_link_libraries(solver ${ARMADILLO_LIBRARIES})
